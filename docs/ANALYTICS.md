# OpenAntigens website analytics

OpenAntigens uses managed Plausible Analytics to measure whether the portal and its downloadable resources are being used. The design uses a small set of aggregate measures that can support maintenance and funder reporting without building visitor profiles.

## Account setup

The Plausible site must be configured for `openantigens.org`. The assigned script is:

```html
<script async src="https://plausible.io/js/pa-Mda6DF7h4_8b17ZUlD_9E.js"></script>
```

The renderer adds the script and initialization code to every page. The script identifier is site-specific but is not a password or API credential. Do not add a Plausible API key, account password, or dashboard-sharing link to the repository.

In Plausible, keep the dashboard private and create these five custom event goals with exact, case-sensitive names:

- `Download PDB Full`
- `Download PDB Selection`
- `Download PNG`
- `Copy TSV`
- `Copy FASTA`

Plausible creates the automatic `File Download` goal after the first tracked file click. The portal limits that automatic measurement to `.tsv` and `.json` files. Outbound-link and form-submission measurement are disabled.

## Event dictionary

| Measure | Trigger | Data attached by OpenAntigens | Interpretation |
| --- | --- | --- | --- |
| Pageview | A tracked page loads | No custom properties | Approximate page traffic |
| `File Download` | A visitor clicks a same-domain `.tsv` or `.json` link | File URL handled by Plausible | A file link was clicked from a tracked page |
| `Download PDB Full` | The full AlphaFold PDB browser download starts | No properties | A full-structure export was initiated |
| `Download PDB Selection` | A selected-region PDB browser download starts | No properties | A selected-region structure export was initiated |
| `Download PNG` | The 3D viewer creates a PNG and starts its download | No properties | A viewer screenshot export was initiated |
| `Copy TSV` | The browser confirms the selected-region TSV was written to the clipboard | No properties | A TSV copy completed |
| `Copy FASTA` | The browser confirms the selected-region FASTA was written to the clipboard | No properties | A FASTA copy completed |

The five custom event names are enforced by an allowlist in the renderer. Unknown names throw an error during development instead of silently expanding collection.

## Deliberate exclusions

OpenAntigens does not attach search terms, protein sequences, residue selections, construct boundaries, clipboard contents, names, email addresses, or other free text to analytics events. The integration does not use analytics cookies, advertising identifiers, session replay, cross-site tracking, form tracking, or an ad-blocker bypass.

Plausible still processes the fields needed for its aggregate dashboard, including page paths, referrers, approximate location, and device, browser, and operating-system categories. See the public `privacy.html` page and [Plausible's data policy](https://plausible.io/data-policy) for the visitor-facing description.

## Production verification

Plausible ignores localhost, loopback addresses, and `file://` pages. Verify analytics only after deploying a build that contains the integration.

1. In Plausible, confirm the site domain is exactly `openantigens.org` and create the five custom event goals above.
2. Open the production site with browser blocking disabled for this test.
3. Load the index, a support page, and a target report.
4. Click one static TSV or JSON download.
5. On a target with a structure, test each PDB/PNG export and each copy action.
6. Confirm the pageviews, `File Download`, and five custom events appear in Plausible's realtime view.
7. Confirm failed screenshot or clipboard actions do not produce their success events.

The static download goal counts clicks from tracked pages. Direct file URLs opened from email, a citation, or another site are not counted. Ad blockers and network failures also reduce measured traffic.

## Quarterly funder report

Use a fixed calendar-quarter window and export or record:

- unique visitors, visits, and pageviews;
- the ten most-viewed pages, separating the index, support pages, and target reports;
- total and unique `File Download` clicks, with the most-clicked file URLs;
- counts for each of the five custom export events;
- top referring sites and country-level reach;
- comparison with the previous quarter, when the previous quarter used the same measurement design.

Label the numbers as approximate. A click or custom export event does not prove that a transfer completed, that a scientist used the file, or that the resource affected an experiment. Keep that distinction in funder materials.

Do not share the live dashboard or visitor-level data with funders. Suppress any geographic, referral, page, or event breakdown with fewer than five unique visitors. Store only the reviewed quarterly aggregate summary used for historical reporting.

## Retention and change control

Plausible analytics remain in the provider account while OpenAntigens uses the service. If OpenAntigens is removed from Plausible, the site administrator must delete its provider-held analytics. Reviewed quarterly aggregate summaries may be retained for historical comparisons.

Any new event, file extension, custom property, or enhanced-measurement option changes the collection boundary. Update the allowlist, tests, this document, and `privacy.html` in the same commit. Review the public wording before deploying the change.
