# LIVIA proxy (Cloudflare Worker)

LIVIA runs in the browser. A few services it uses cannot be called directly from a web page: MIST (flyrnai.org) sends no
cross-origin (CORS) headers, some hosts of prediction bundles (for example OSF) do the same, and BioGRID needs an access key
that cannot be placed in public page code. `worker.js` is a small proxy for these cases. It:

- answers GET requests only;
- grants cross-origin access only to the origins in `ALLOWED_ORIGINS` (the hosted LIVIA and local copies on ports 8000 and 8765);
- forwards `?url=` requests only to the hosts in `ALLOWED_TARGET_HOSTS`, so it is not an open relay;
- adds the BioGRID key on the server for `?biogrid=1` requests (the key is a Worker secret, never in this file).

Scoring and every residue-level view in LIVIA work without the proxy; only reported-interaction panels and bundles on hosts
without CORS headers need it.

## Asking for a host to be added

Open an issue on https://github.com/flyark/LIVIA with the host name and an example bundle URL. Please do not ask by e-mail;
an issue keeps the allow-list public.

## Running your own copy (independent installations)

1. Fork LIVIA.
2. Install Wrangler (`npm install -g wrangler`) and sign in to your Cloudflare account (`wrangler login`).
3. In `worker.js`, set `ALLOWED_ORIGINS` to the address(es) where your copy of LIVIA is served, and `ALLOWED_TARGET_HOSTS`
   to the hosts you need.
4. Optional, for BioGRID: `wrangler secret put BIOGRID_KEY` (key from https://wiki.thebiogrid.org/doku.php/biogridrest).
5. From this folder: `wrangler deploy`. Wrangler prints the worker's address, `https://<name>.<account>.workers.dev`.
6. Put that address in `LIVIA_PROXY` at the top of `js/livia-core.js` (the only place LIVIA reads it), or define
   `window.LIVIA_CONFIG = { proxy: 'https://<name>.<account>.workers.dev' }` in a script loaded before `js/livia-core.js`.

The free Workers plan covers the `?url=` and `?biogrid=1` routes; the `?fpSummary=` route needs a paid plan (see `wrangler.toml`).
