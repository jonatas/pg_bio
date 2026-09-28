---
name: pg_bio_blog_writer
description: Standardizes the creation of SEO-optimized blog posts for pg_bio discoveries, ensuring they include SQL, math (vector distances), biological learnings, and local 3Dmol.js previews.
---

# pg_bio Discovery Blog Writer

Use this skill whenever documenting a new protein discovery made via `pg_bio` into the `ideia.me` Jekyll blog. 

## Requirements for Every Post

### 1. SEO & Frontmatter
* Create a catchy, highly searchable title (e.g., "Mining Bioplastic Enzymes in the Dark Proteome using AI and SQL").
* Include strong keywords in the Jekyll `categories:` tag (e.g., `bioinformatics`, `pgvector`, `machine-learning`, `structural-biology`).
* Use `<!--more-->` after the introduction to ensure the blog index preview works perfectly.

### 2. The Biological Hook, Ecosystem, & Behavior (The "Why & Where")
* **MANDATORY RESEARCH:** You must look up the exact taxonomy/organism for both the Bait and the Discovery (e.g., using UniProt API or web search) to understand their ecosystem.
* Detail the **ecosystem and behavior**: Where do these organisms live? (e.g., Deep-sea hydrothermal vents, hypersaline lakes, acidic hot springs). How does the organism survive? What is the ecological role of this protein?
* Explain the mechanism of action: How does the protein physically behave to achieve its function?
* Contrast the Bait's environment with the Discovery's environment to highlight the evolutionary divergence.

### 3. The Math & The SQL
* Always include the exact `pg_bio` SQL query used to find the protein. Highlight the vector distance `<=>` and the hybrid sequence operator `<~>`.
* Include a Markdown table showing the Bait, the Discovery, the Vector Distance, and the Hybrid Score.
* Explain the vector distance. (e.g., "A distance of 0.068 means the 3D backbone is almost mathematically identical").

### 4. The 3Dmol.js Interactive Preview
* **CRITICAL RULE:** Do NOT link to external AlphaFold or PDB URLs in the 3Dmol.js viewer. 
* You MUST use curl to download representative `.pdb` files from RCSB PDB (e.g., `https://files.rcsb.org/download/XXXX.pdb`) and save them to `../ideia.me/assets/pdb/`.
* Configure the 3D viewer in the HTML block to load the local assets (`data-href="/assets/pdb/your_file.pdb"`).
* Put the Bait and the Discovery side-by-side using the flexbox layout, colored differently (e.g., cyan vs magenta), with a caption pointing out an interesting structural feature (like a binding pocket).

### 5. The Learning/Conclusion
* Synthesize what this means for computational biology. Remind the reader that this was done in seconds using PostgreSQL, bypassing months of wet-lab work.

### 6. Interactive 3Dmol.js Features (Built-in)
* Mention to the user that the blog features an interactive 3D plugin (`pg_bio_sync.js`).
* The user can double-click either 3D viewer to lock their cameras together (synchronized rotation/tilt).
* The user can click any fragment on one protein, and it will automatically highlight the matching residue (across all complex chains!) on the opposite protein by swapping colors.
