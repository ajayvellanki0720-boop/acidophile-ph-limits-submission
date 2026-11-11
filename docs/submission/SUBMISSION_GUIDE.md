# PeerJ upload package

Upload the final `Acidophile_PeerJ_Line_Numbered_Manuscript.docx` as the manuscript. It contains native continuous Word line numbering, a structured abstract and figure legends. Figures are separate, as required by the journal upload instructions.

Upload `figures/Figure_1.png` through `Figure_4.png` as the four figures. These 600-dpi files exceed the portal's 900-pixel minimum in both dimensions. SVG sources are retained in the repository; the portal lists PDF/EPS for vector uploads. Figure 4 combines panels A and B. Upload only one format for each figure unless requested otherwise.

The portal rejected the valid Word file with an unexpected request for an EPS extension. Its accepted RTF alternative, `Acidophile_PeerJ_Line_Numbered_Manuscript.rtf`, was exported with LibreOffice and uploaded successfully. The Word source remains the editable master. The current master has a separate author cover page and no headers or footers, following the portal's manuscript checks. The portal requires a short title even though the manuscript has no running title.

To regenerate the RTF with LibreOffice installed, run `soffice --headless --convert-to rtf --outdir docs/submission docs/submission/Acidophile_PeerJ_Line_Numbered_Manuscript.docx` after rebuilding the Word source.

Upload `Supplemental_Data_S1.zip` as Supplemental Data S1, titled “Curated inputs and stored analysis outputs for extracytoplasmic surface charge and reported minimum growth pH.” It is about 21 MB, below the stated 30 MB individual-file limit. This is intended as journal-hosted supplemental data; the GitHub repository supplies code, figures, retrieval instructions and a public project record. The code and reproducibility archive is preserved at https://doi.org/10.5281/zenodo.24147308 (version 1.0.0, commit 97257f22dde358b920c70a67aa921a08556b41fb). Include this DOI in the submission system data availability declaration.

Repository: https://github.com/ajayvellanki0720-boop/acidophile-ph-limits-submission

## Author declaration confirmed

The AI disclosure records Anthropic Claude Opus 5 (version supplied by the author) and OpenAI Codex using GPT-6.1, their uses and the distinction between data-driven plotting and generated imagery. On 5 November 2026 the author confirmed the following statement, which is reflected in the manuscript:

“I have reviewed the manuscript, references, code and data visualisations and confirm the originality and accuracy of the content. I have checked the applicable terms of the AI tools used and confirm their suitability for this publication. I have the right to publish the figures and other included material, have obtained any required permissions, and take responsibility for the integrity of the complete article, including the references.”

Retain original and revised drafts, relevant prompts and figure versions for editorial requests. The available package does not contain a complete export of earlier AI conversations; those records must be retained by the author where available.

## Scientific scope of this revision

The revised wording reports observational associations with reported growth minima. It discloses the neighbour-joining tree method and the absence of a full-sample maximum-likelihood sensitivity analysis. No new phylogeny or biological experiment was performed during this revision. The canonical panel remains distinct from its deterministic rebuild, with the measured headline sensitivity stated explicitly. Historical documents preserve earlier analyses and language; the final submission document is the current manuscript.

Current instructions consulted: https://peerj.com/journals/peerj/author-instructions/ and https://peerj.com/journals/peerj/policies/#generative-ai
