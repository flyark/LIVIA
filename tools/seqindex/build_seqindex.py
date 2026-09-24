#!/usr/bin/env python3
"""
build_seqindex.py - build the static sequence index that chain identification looks up first.

Plain-text shards served by GitHub Pages; the browser fetches one small file per lookup.

  x/<hhh>.txt   Exact index. Key = SHA-256 (hex) of the uppercase sequence; <hhh> = its first
                3 hex characters (4,096 files). Row: key[3:16]  accession  gene  taxid  s|t  length  afdb
                (s = Swiss-Prot reviewed, t = unreviewed; isoforms keep their -N suffix; afdb = 1 when
                UniProt cross-references an AlphaFold DB model for the accession, 0 when it does not).
                Holds Swiss-Prot + isoforms, the model organisms' reference proteomes, and their
                unreviewed entries outside the reference proteome that have an AFDB model (a gene's
                reference entry can lack a model that an identical-sequence sibling has).
  s/<hhhh>.txt  Seed index, for fragments and point mutants (model organisms, one sequence per
                gene). A 12-mer is a seed when the low 4 bits of H1 are zero, so a query picks the
                same seeds as its parent. File = bits 4-17 of H1 (16,384 files). Row: H2  pid  pos (hex).
  p/<hhh>.txt   Seed proteins, 256 per file: file = pid >> 8, line = pid & 255.
                Row: accession  gene  taxid  s|t  length  afdb
  t/<nn>.txt    Organism names, file = taxid % 100. Row: taxid  scientific name
  meta.json     UniProt release, counts, format version.

H1 and H2 are 32-bit FNV-1a over the 12 ASCII letters (two offset bases) finished with the
murmur3 fmix32 avalanche: cheap, and bit-identical in JavaScript (Math.imul) and numpy (uint32).

Usage: build_seqindex.py --cache DIR --out DIR     (downloads the UniProt files into --cache)
"""
import argparse, collections, datetime, gzip, hashlib, json, pathlib, re, sys, urllib.parse, urllib.request
import numpy as np

UNIPROT = 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase'
REST = 'https://rest.uniprot.org/uniprotkb/stream?compressed=true&'
ORGANISMS = {   # taxid: (reference proteome, kingdom folder)
    '9606': ('UP000005640', 'Eukaryota'), '10090': ('UP000000589', 'Eukaryota'), '10116': ('UP000002494', 'Eukaryota'),
    '7227': ('UP000000803', 'Eukaryota'), '559292': ('UP000002311', 'Eukaryota'), '284812': ('UP000002485', 'Eukaryota'),
    '6239': ('UP000001940', 'Eukaryota'), '7955': ('UP000000437', 'Eukaryota'), '3702': ('UP000006548', 'Eukaryota'),
    '83333': ('UP000000625', 'Bacteria'),
}
K, SAMPLE_BITS, SEED_SHARD_BITS, MAX_OCC = 12, 4, 14, 60
FNV_PRIME, BASIS1, BASIS2 = 16777619, 2166136261, 2166136261 ^ 0x5BD1E995
HDR = re.compile(r'^>(sp|tr)\|([^|]+)\|\S+\s.*?\sOS=(.*?)\sOX=(\d+)(?:.*?GN=(\S+))?')


def fetch(url, dest):
    if not dest.exists():
        print('  downloading', url, file=sys.stderr)
        tmp = dest.with_suffix('.part')
        urllib.request.urlretrieve(url, tmp)
        tmp.rename(dest)
    return dest


def read_fasta(path):
    head, buf = None, []
    with gzip.open(path, 'rt') as f:
        for line in f:
            if line.startswith('>'):
                if head:
                    yield head, ''.join(buf)
                head, buf = line.rstrip(), []
            else:
                buf.append(line.strip())
    if head:
        yield head, ''.join(buf)


def records(path):
    for head, seq in read_fasta(path):
        m = HDR.match(head)
        if m:
            yield dict(db='s' if m.group(1) == 'sp' else 't', acc=m.group(2), org=m.group(3), tax=m.group(4),
                       gene=m.group(5) or '', seq=seq.upper())


def fmix32(h):
    h ^= h >> np.uint32(16); h *= np.uint32(0x85EBCA6B)
    h ^= h >> np.uint32(13); h *= np.uint32(0xC2B2AE35)
    h ^= h >> np.uint32(16)
    return h


def fnv_windows(arr, starts, basis):
    """FNV-1a + fmix32 of the K-mer starting at each index in `starts` of the uint8 array `arr`."""
    h = np.full(starts.shape, basis, dtype=np.uint32)
    for j in range(K):
        h ^= arr[starts + j].astype(np.uint32)
        h *= np.uint32(FNV_PRIME)
    return fmix32(h)


def afdb_flag(acc, afdb):
    return '1' if acc.split('-')[0] in afdb else '0'     # an isoform row carries its entry's flag


def build_exact(recs, out, taxnames, afdb):
    shards = collections.defaultdict(set)
    for r in recs:
        key = hashlib.sha256(r['seq'].encode()).hexdigest()
        shards[key[:3]].add(f"{key[3:16]}\t{r['acc']}\t{r['gene']}\t{r['tax']}\t{r['db']}\t{len(r['seq'])}\t{afdb_flag(r['acc'], afdb)}")
        taxnames.setdefault(r['tax'], r['org'])
    d = out / 'x'; d.mkdir(parents=True, exist_ok=True)
    rows = 0
    for i in range(4096):
        name = f'{i:03x}'
        lines = sorted(shards.get(name, ()))
        rows += len(lines)
        (d / f'{name}.txt').write_text(''.join(l + '\n' for l in lines))
    return rows


def build_seeds(prots, out, afdb):
    d = out / 's'; d.mkdir(parents=True, exist_ok=True)
    pd = out / 'p'; pd.mkdir(parents=True, exist_ok=True)
    for base in range(0, len(prots), 256):
        chunk = prots[base:base + 256]
        (pd / f'{base >> 8:03x}.txt').write_text(''.join(
            f"{r['acc']}\t{r['gene']}\t{r['tax']}\t{r['db']}\t{len(r['seq'])}\t{afdb_flag(r['acc'], afdb)}\n" for r in chunk))
    shard_l, key_l, pid_l, pos_l = [], [], [], []
    CHUNK = 4_000_000
    i = 0
    while i < len(prots):
        group, total = [], 0
        while i < len(prots) and (total < CHUNK or not group):
            group.append(i); total += len(prots[i]['seq']) + 1; i += 1
        text = b'\0'.join(prots[p]['seq'].encode() for p in group) + b'\0'
        arr = np.frombuffer(text, dtype=np.uint8)
        owner = np.repeat(np.array(group, dtype=np.int64), [len(prots[p]['seq']) + 1 for p in group])
        offsets = np.concatenate([[0], np.cumsum([len(prots[p]['seq']) + 1 for p in group])[:-1]])
        pos_in = np.arange(len(arr), dtype=np.int64) - np.repeat(offsets, [len(prots[p]['seq']) + 1 for p in group])
        n = len(arr) - K + 1
        starts = np.arange(n, dtype=np.int64)
        seps = np.concatenate([[0], np.cumsum(arr == 0)])
        valid = (seps[starts + K] - seps[starts]) == 0
        starts = starts[valid]
        h1 = fnv_windows(arr, starts, BASIS1)
        keep = (h1 & np.uint32((1 << SAMPLE_BITS) - 1)) == 0
        starts, h1 = starts[keep], h1[keep]
        h2 = fnv_windows(arr, starts, BASIS2)
        shard_l.append((h1 >> np.uint32(SAMPLE_BITS)) & np.uint32((1 << SEED_SHARD_BITS) - 1))
        key_l.append(h2); pid_l.append(owner[starts]); pos_l.append(pos_in[starts])
    shard, key = np.concatenate(shard_l), np.concatenate(key_l)
    pid, pos = np.concatenate(pid_l), np.concatenate(pos_l)
    order = np.lexsort((pos, pid, key, shard))
    shard, key, pid, pos = shard[order], key[order], pid[order], pos[order]
    # drop over-represented seeds (low complexity / repeats): useless for placement and bulky
    grp = np.concatenate([[True], (shard[1:] != shard[:-1]) | (key[1:] != key[:-1])])
    gid = np.cumsum(grp) - 1
    counts = np.bincount(gid)
    ok = counts[gid] <= MAX_OCC
    shard, key, pid, pos = shard[ok], key[ok], pid[ok], pos[ok]
    bounds = np.searchsorted(shard, np.arange((1 << SEED_SHARD_BITS) + 1))
    for s in range(1 << SEED_SHARD_BITS):
        a, b = bounds[s], bounds[s + 1]
        (d / f'{s:04x}.txt').write_text(''.join(
            f'{k:08x}\t{p:x}\t{q:x}\n' for k, p, q in zip(key[a:b].tolist(), pid[a:b].tolist(), pos[a:b].tolist())))
    return len(key)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cache', required=True, type=pathlib.Path)
    ap.add_argument('--out', required=True, type=pathlib.Path)
    a = ap.parse_args()
    a.cache.mkdir(parents=True, exist_ok=True); a.out.mkdir(parents=True, exist_ok=True)
    sprot = fetch(f'{UNIPROT}/complete/uniprot_sprot.fasta.gz', a.cache / 'uniprot_sprot.fasta.gz')
    iso = fetch(f'{UNIPROT}/complete/uniprot_sprot_varsplic.fasta.gz', a.cache / 'uniprot_sprot_varsplic.fasta.gz')
    rel = urllib.request.urlopen(f'{UNIPROT}/complete/reldate.txt').read().decode()
    release = (re.search(r'Release (\d{4}_\d{2})', rel) or [None, 'unknown'])[1]
    rest = lambda query, fmt: REST + urllib.parse.urlencode({'query': query, 'format': fmt})
    lists = [fetch(rest('(reviewed:true) AND (database:alphafolddb)', 'list'), a.cache / 'afdb_sprot.list.gz')]
    per_gene, additional, outside = {}, {}, {}
    for tax, (up, kingdom) in ORGANISMS.items():
        base = f'{UNIPROT}/reference_proteomes/{kingdom}/{up}/{up}_{tax}'
        per_gene[tax] = fetch(base + '.fasta.gz', a.cache / f'{up}_{tax}.fasta.gz')
        additional[tax] = fetch(base + '_additional.fasta.gz', a.cache / f'{up}_{tax}_additional.fasta.gz')
        lists.append(fetch(rest(f'(proteome:{up}) AND (database:alphafolddb)', 'list'), a.cache / f'afdb_{up}.list.gz'))
        outside[tax] = fetch(rest(f'organism_id:{tax} AND database:alphafolddb AND NOT proteome:{up} AND reviewed:false', 'fasta'),
                             a.cache / f'outside_{up}_{tax}.fasta.gz')
    afdb = set()
    for p in lists:
        with gzip.open(p, 'rt') as f:
            afdb.update(line.strip() for line in f if line.strip())
    for p in outside.values():                          # selected by having a model
        afdb.update(r['acc'] for r in records(p))
    print(f'AFDB-backed accessions: {len(afdb):,}', file=sys.stderr)

    taxnames = {}
    exact_src = [sprot, iso] + list(per_gene.values()) + list(additional.values()) + list(outside.values())
    exact_rows = build_exact((r for p in exact_src for r in records(p)), a.out, taxnames, afdb)
    print(f'exact index: {exact_rows:,} rows', file=sys.stderr)

    seen, prots = set(), []
    for tax in ORGANISMS:
        for r in records(per_gene[tax]):
            if r['seq'] not in seen and len(r['seq']) >= K:
                seen.add(r['seq']); prots.append(r)
    seeds = build_seeds(prots, a.out, afdb)
    print(f'seed index: {len(prots):,} proteins, {seeds:,} seeds', file=sys.stderr)

    td = a.out / 't'; td.mkdir(exist_ok=True)
    buckets = collections.defaultdict(list)
    for tax, name in taxnames.items():
        buckets[int(tax) % 100].append(f'{tax}\t{name}')
    for b in range(100):
        (td / f'{b:02d}.txt').write_text(''.join(l + '\n' for l in sorted(buckets.get(b, []))))
    (a.out / 'meta.json').write_text(json.dumps({
        'format': 1, 'uniprot_release': release, 'built': datetime.date.today().isoformat(),
        'exact': {'rows': exact_rows, 'files': 4096, 'key': 'sha256 hex, first 3 = file, next 13 = row key'},
        'seed': {'proteins': len(prots), 'seeds': int(seeds), 'k': K, 'sample': 1 << SAMPLE_BITS,
                 'files': 1 << SEED_SHARD_BITS, 'max_occurrences': MAX_OCC},
        'organisms': sorted(ORGANISMS, key=int),
    }, indent=1) + '\n')


if __name__ == '__main__':
    main()
