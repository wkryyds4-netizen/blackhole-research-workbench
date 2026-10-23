# Contributing

Please separate numerical methods, dataset adapters and UI changes. Include
units, reference frames, source versions and scientific limitations with new
data. Do not label archived data as held-out validation without verifying model
training membership.

Run `npm ci`, then `npm test` for the portable frontend contracts. For browser
checks, install Playwright Chromium with `npx playwright install chromium`,
serve `simulator/public` on loopback port 8765, set
`PLAYWRIGHT_CHANNEL=chromium`, and run `npm run test:release`.

Large raw datasets, Python environments, downloaded third-party repositories,
credentials and data lacking redistribution permission must stay out of Git.
Use fixed source URLs and checksums instead. The public dataset catalog is
`simulator/public/datasets/catalog.public.json`; the ignored local `catalog.json`
may additionally refer to privately held research snapshots.

Numerical changes need tests against independent calculations or known
synthetic cases. A lower mismatch alone does not establish improved physics.
