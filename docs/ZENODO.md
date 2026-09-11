# Mint a Zenodo DOI from mineralhoiland/QVC

Zenodo’s GitHub integration **cannot see a private repository**. Make `QVC`
public first, connect Zenodo, then mint the DOI from a GitHub Release.

Do this from the snapshot git, not the QVCCursor lab tree:

```bash
cd /Users/mineralhoiland/Code/QVCCursor/releases/qvc
```

## 1. Push the current snapshot

Include `.zenodo.json`, `CITATION.cff`, and the revised paper sources.

```bash
git status
git add -A
git commit -m "Prepare v0.4.1 for Zenodo: metadata, ORCID, public manuscript."
git push origin main
```

## 2. Make the GitHub repo public

On GitHub: **Settings → General → Danger zone → Change repository visibility → Public**.

Or:

```bash
gh repo edit mineralhoiland/QVC --visibility public --accept-visibility-change-consequences
```

Wait until https://github.com/mineralhoiland/QVC loads while logged out.

## 3. Connect Zenodo to GitHub

1. Open https://zenodo.org and **Log in with GitHub** (same account as `mineralhoiland`).
2. Open https://zenodo.org/account/settings/github/
3. Click **Sync now** if `QVC` is missing.
4. Flip the switch **on** for `mineralhoiland/QVC`.

Zenodo only lists **public** repos. It will install a webhook.

## 4. Create the GitHub Release (this mints the DOI)

Do **not** tag until the Zenodo switch is on. A release made while the hook is off is not archived.

Do **not** retag `v0.4.0` or `release`. Zenodo already marked those as received (HTTP 409) and minted no public DOI. After metadata is on `main`, use a **new** tag:

```bash
git tag -a v0.4.1 -m "QVC preprint compilation v0.4.1"
git push origin v0.4.1
gh release create v0.4.1 \
  --title "v0.4.1 — QVC preprint compilation" \
  --notes-file docs/RELEASE_NOTES.md
```

Wait one or two minutes, then open the GitHub tab on Zenodo. You should see:

- a **version DOI** for `v0.4.1` (cite this snapshot)
- a **concept DOI** that always points at the latest version

## 5. Put the DOI back into the paper and CITATION.cff

1. Copy both DOIs from the Zenodo record.
2. Paste the version DOI into `CITATION.cff` (`doi:`) and into Appendix “Data and Code Availability” in `paper/QVC_arXiv_v4.tex` (and the lab copy `tex/QVC_arXiv_v4/`).
3. Add a README badge using the **concept DOI** if you want it to track later versions.
4. Commit the DOI onto `main`. Tag a later patch only if you need the badge inside the archived tree.

Example badge:

```markdown
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
```

## Dual license

Zenodo’s GitHub deposit uses one license field. `.zenodo.json` sets **CC BY 4.0** for the compilation as a whole. The notes field states that `compute/` and `sim/` remain MIT. Do not change that after the DOI is minted unless you issue a new version.

## If you must keep GitHub private

The GitHub hook will not work. Zip `releases/qvc` (without `.git` and without `node_modules`), upload at https://zenodo.org/uploads/new as a preprint, paste the same metadata, and publish. Prefer the GitHub path; that is how TABD was minted.

## What this does not do

A Zenodo DOI is not an arXiv identifier and is not peer review. After the DOI exists, submit the TeX bundle to arXiv and cite the version DOI in the comments line and in Appendix `app:data`.
