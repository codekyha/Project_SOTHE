# Release guide: GitHub branch `p2`, tag `p2-v1.0.0`, and an own Zenodo record

This guide takes SOTHE-P2 from this folder to three places:

- a GitHub branch `p2` of `codekyha/Project_SOTHE`, holding exactly this folder;
- a tag `p2-v1.0.0`, which the manuscript cites as `[RELEASE-TAG]`;
- a Zenodo record of its own (manual upload, not the GitHub integration), whose version DOI fills `[ZENODO-VERSION-DOI]`.

Every git, GitHub and Zenodo step below is done by the PI (HR-10). The tools in `tools/` only build files locally.

**Order matters.** Reserve the DOI first, write it into the files, then commit and tag, then upload the tarball of the same tree. The tag, the Zenodo files and the DOI then all describe the same bytes.

## 0. Before releasing

```bash
python -m unittest discover -s tests          # all 41 tests pass (about 45 s)
python p2.py run      # or, on Altay: python p2.py submit   (a full run; keep its results tarball for Zenodo)
python p2.py show     # SUMMARY.md
```

The full run is expected to show:

| Stage | Expected |
|---|---|
| S1 | G4a all pass |
| S2 | 9/9 converged, scaling ratio within 5 %, recoil deviation within 15 % |
| S3 | run A reproduces the archive |
| S4 | 29/29 and 35/35 |
| S5 | five figures (PDF, EPS, PNG) and their captions |
| S6 | 29/29, ten CSV files, six figures with captions, parity with MATLAB at round-off |
| S7 | 300 configurations |

Four files carry the version and must agree: `VERSION`, `sothe_p2/__init__.py` (`__version__`), `pyproject.toml` (`version`), and `CITATION.cff` (`version`, `date-released`). `CHANGELOG.md` needs an entry for the new version. `tools/build_release.py` refuses to build when the four files disagree.

## 1. Build the release tarball

```bash
python tools/build_release.py            # MANIFEST.sha256, dist/SOTHE-P2_v1.0.0.tar.gz and .sha256
python tools/build_release.py --check    # MANIFEST against the files: 0 changed, 0 unlisted
```

The tarball is reproducible: the same tree always gives the same bytes. It unpacks into one folder, `SOTHE-P2/`. It is also the file to upload to Altay ([README_UHeM_Altay.md](README_UHeM_Altay.md)).

## 2. Reserve the Zenodo version DOI

1. <https://zenodo.org>: log in with the account that owns the Paper-1 record (doi:10.5281/zenodo.20713660).
2. **New upload.** Do not open the Paper-1 record and choose "New version": that would put SOTHE-P2 under Paper 1's concept DOI.
3. Under **Digital Object Identifier**, answer "Do you already have a DOI?" with *No*, then press **Get a DOI now!**. Copy the DOI, for example `10.5281/zenodo.1234567`.
4. **Save draft**. Do not publish yet.
5. Write the DOI into the files and rebuild:

   ```bash
   python tools/build_release.py --doi 10.5281/zenodo.1234567
   ```

   This fills the `doi:` line of `CITATION.cff` and the "Zenodo DOI of this version" line of `README.md`, then rebuilds `MANIFEST.sha256` and the tarball.

## 3. GitHub: branch `p2`, commit, tag

`p2` is a standalone branch. An orphan branch keeps the history of `main` (Paper 1) out of it. Use Git Bash on Windows, or any shell on Linux or macOS, in a fresh folder:

```bash
git clone https://github.com/codekyha/Project_SOTHE.git Project_SOTHE-p2
cd Project_SOTHE-p2
git switch --orphan p2                   # Git >= 2.23; older: git checkout --orphan p2 && git rm -rf .
git status                               # nothing tracked; delete any untracked leftovers
tar -xzf /path/to/dist/SOTHE-P2_v1.0.0.tar.gz --strip-components=1
git add -A
git status                               # the pack files only: runs/, dist/, config/site.env are in .gitignore
git commit -m "SOTHE-P2 1.0.0: Python rewrite of the Paper-2 computations (P2_altay_pack 1.0.0, 68 MATLAB files)"
git push -u origin p2
git tag -a p2-v1.0.0 -m "SOTHE-P2 1.0.0"
git push origin p2-v1.0.0
```

PowerShell uses the same commands; `tar` is built into Windows 10 and later. `.gitattributes` makes Git check out LF line endings, so `python p2.py verify` also passes in a Windows clone.

To check what was pushed:

```bash
git clone -b p2-v1.0.0 https://github.com/codekyha/Project_SOTHE.git check-p2 && cd check-p2 && python p2.py verify
```

### 3a. GitHub release (optional)

**Check the Zenodo-GitHub integration first.** On zenodo.org, open the user menu and choose GitHub. If `codekyha/Project_SOTHE` is switched *on* there, publishing a GitHub *release* makes Zenodo archive it automatically, as a second record with another DOI. For this release, either switch the repository off first, or skip the GitHub release: the tag alone is what the manuscript cites.

On the web:

1. Go to Releases, then Draft a new release.
2. Tag `p2-v1.0.0`, target `p2`, title `SOTHE-P2 1.0.0`.
3. Notes: the 1.0.0 section of `CHANGELOG.md`.
4. Attach `dist/SOTHE-P2_v1.0.0.tar.gz` and its `.sha256`, then Publish.

With the GitHub CLI:

```bash
gh release create p2-v1.0.0 dist/SOTHE-P2_v1.0.0.tar.gz dist/SOTHE-P2_v1.0.0.tar.gz.sha256 \
   --repo codekyha/Project_SOTHE --target p2 --title "SOTHE-P2 1.0.0" --notes "See CHANGELOG.md (1.0.0)."
```

GitHub also offers its own "Source code" archives. They hold the same files, but in another top folder and without `dist/`. The attached `SOTHE-P2_v1.0.0.tar.gz` is the one whose checksum is published.

## 4. The Zenodo pack

Build the pack from the committed tree (`git status` clean), with the reserved DOI and the results of the run that supports the paper:

```bash
python tools/make_zenodo_pack.py --doi 10.5281/zenodo.1234567 \
       --results /path/to/SOTHE-P2_results_<RUNID>.tar.gz     # repeat --results for more runs; --extra for other files
```

It writes `dist/zenodo_SOTHE-P2_v1.0.0/`:

| File | Content |
|---|---|
| `SOTHE-P2_v1.0.0.tar.gz` (+ `.sha256`) | the release, identical to the tag |
| `SOTHE-P2_results_<RUNID>.tar.gz` (+ `.sha256`) | each `--results` (a run folder is packed reproducibly) |
| `README_ZENODO.md` | description of the record and its files |
| `zenodo_metadata.json` | every field of the form (also valid for the Zenodo REST API) |
| `CHECKSUMS.sha256` | SHA-256 of every file |
| `UPLOAD_STEPS.txt` | the steps below |

Which run goes into the record is the PI's decision. Under HR-3, a Python run is RECORDED evidence; a canonical MATLAB R2025b output archive can be added with `--extra`.

## 5. Upload and publish on Zenodo

In the draft saved at step 2:

1. **Files.** Drag in every file of `dist/zenodo_SOTHE-P2_v1.0.0/` except `UPLOAD_STEPS.txt`.
2. **Basic information**, copied from `zenodo_metadata.json`:
   - resource type *Software*;
   - title;
   - publication date;
   - creator *Oguz, Hasan*, with ORCID 0000-0001-7484-4415 and affiliation Istanbul Okan University; Pamukkale University;
   - description: paste it, or use the text of `README_ZENODO.md`;
   - version `1.0.0`;
   - licence MIT;
   - keywords.
3. **Related works**:
   - the tag URL `https://github.com/codekyha/Project_SOTHE/tree/p2-v1.0.0`, as *Is supplement to*;
   - `10.5281/zenodo.20713660` (Ref. [25] code and data), as *References*;
   - `10.1088/1361-6382/ae811b` (Ref. [25] article), as *References*.
4. **Preview**, then **Publish**.

The record then has two DOIs:

- the **version DOI**, the one reserved at step 2, which cites exactly these files;
- the **concept DOI**, which always resolves to the newest version.

The manuscript cites the version DOI.

A different title is possible: the Paper-1 record is titled "Supplementary: code and data for '<paper title>'". The same pattern can be used once the Paper-2 title is final.

## 6. After publication (PI, B3)

- **Manuscript.**
  - `[RELEASE-TAG]` becomes `p2-v1.0.0`.
  - `[ZENODO-VERSION-DOI]` becomes the version DOI, in the form the data-availability text uses (for example `\href{https://doi.org/10.5281/zenodo.1234567}{doi:10.5281/zenodo.1234567}`).
  - Record the fill in the DELTA log.
- **Repository.** Nothing more to change. The DOI is already in `CITATION.cff` and `README.md` (step 2).
- **When Paper 2 is published.** Edit the Zenodo record, add the article DOI under Related works (*Is supplement to*), and publish. A metadata edit keeps the DOI.

## 7. Next versions

1. Update `VERSION`, `sothe_p2/__init__.py`, `pyproject.toml`, `CITATION.cff` (`version`, `date-released`) and `CHANGELOG.md`.
2. On Zenodo, open the SOTHE-P2 record, choose **New version**, and reserve the new version DOI.
3. Run `python tools/build_release.py --doi <new DOI>` (it replaces the previous DOI in `CITATION.cff` and `README.md`), commit on `p2`, tag `p2-v1.0.1`, and push the branch and the tag.
4. Run `python tools/make_zenodo_pack.py --doi <new DOI> ...`, replace the files in the new-version draft, and Publish.

The concept DOI stays the same; each version keeps its own DOI.
