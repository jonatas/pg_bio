---
name: pg_bio_blog_writer
description: Standardizes the creation of highly enthusiastic, SEO-optimized blog posts for pg_bio discoveries, ensuring they include SQL, math (vector distances), rich biological context, practical applications, and local 3Dmol.js previews.
---

# pg_bio Discovery Blog Writer

Use this skill whenever documenting a new protein discovery made via `pg_bio` into the `ideia.me` Jekyll blog. 

## Requirements for Every Post

### 0. NO STATIC TEMPLATES (Dynamic Reasoning Required)
* **CRITICAL:** Do NOT use static text templates (e.g., "Unearthing [Family]..." or "As the pg_bio autonomous night pipeline continues...").
* You must use your reasoning capabilities to analyze the specific findings (the proteins, the organisms, the vector distances) and dynamically craft original, flowing prose for every single post. 
* Write naturally in a blog style that sounds like a real human (the author) sharing an exciting new discovery, adapting the structure and narrative to fit the actual biology of the discovery rather than filling in blanks in a pre-written paragraph.

### 1. Tone, SEO & Frontmatter
* **Enthusiastic Tone:** The writing must be passionate and enthusiastic! Highlight the excitement of using AI and SQL to bypass months of wet-lab work and uncover hidden secrets of the dark proteome. Make it sound organic and tailored to the specific discovery.
* Create a catchy, highly searchable title that is unique to the specific finding.
* Include strong keywords in the Jekyll `categories:` tag (e.g., `bioinformatics`, `pgvector`, `machine-learning`, `structural-biology`).
* Use `<!--more-->` after the introduction to ensure the blog index preview works perfectly.

### 2. Deep Dive: Dynamic Reasoning on Biology
* **Analyze the Bait:** Do not just list the function; weave a narrative around *why* the bait's role is fascinating in its native ecosystem. What makes it biologically significant?
* **Analyze the Discovery (Orphan):** Use reasoning to bridge the gap between the well-studied bait and the unknown orphan. Contrast their environments, taxonomy, or predicted structures.
* **MANDATORY RESEARCH & INSIGHT:** You must look up the exact taxonomy/organism and function for both the Bait and the Discovery (e.g., using UniProt API or web search). Do not just state facts—synthesize them into a compelling hypothesis about how the orphan might have evolved or adapted.

### 3. Practical Applications: Extrapolating the Future
* Instead of generic suggestions (like a bulleted list of biotechnology/medicine), dynamically reason about *specific* real-world applications based on the organism's extremophile traits or unique environment.
* Tell a story about how discovering this exact protein variant could solve a real problem (e.g., a specific industrial hurdle, a new bioremediation technique), linking it back to the unique properties implied by its origin.

### 4. The Math & The SQL
* Always include the exact `pg_bio` SQL query used to find the protein. Highlight the vector distance `<=>` and the hybrid sequence operator `<~>`.
* Include a Markdown table showing the Bait, the Discovery, the Vector Distance, and the Hybrid Score.
* Explain the vector distance. (e.g., "A distance of 0.068 means the 3D backbone is almost mathematically identical!").

### 5. Interactive 3Dmol.js Validation & Chemistry Teaching
* **CRITICAL RULE:** Do NOT link to external AlphaFold or PDB URLs in the 3Dmol.js viewer. 
* You MUST use curl to download representative `.pdb` files from RCSB PDB or AF DB and save them to `/Users/jonatas/code/ideia.me/assets/models/`.
* Configure the 3D viewer in the HTML block to load the local assets (`data-href="/assets/models/your_file.pdb"`).
* Put the Bait and the Discovery side-by-side using the flexbox layout, colored differently (e.g., cyan vs magenta).
* **MANDATORY TEACHING INJECTION:** You MUST include an interactive `<script>` block and HTML UI below the viewers. This UI must have buttons (e.g., "Highlight Conserved Core", "Highlight Adaptations") that execute JS to manipulate the 3Dmol.js viewer (`$3Dmol.viewers`).
* The JS must use methods like `.setStyle()` and `.addLabel()` to highlight specific chemical features (e.g., active sites, hydrophobic cores, acidic surfaces).
* **DYNAMIC LEARNING PARAGRAPHS:** When a user clicks a button, a hidden text `div` below the buttons must dynamically update with a paragraph explaining *what the chemistry means* (e.g., explaining what a conserved active site is, or why acidic shells help in salt flats). Always use this to teach real chemistry!

### 6. The Learning/Conclusion
* Synthesize what this means for computational biology. Remind the reader of the immense power of native PostgreSQL multiomics engines scanning millions of vectors in milliseconds.

### 7. Interactive 3Dmol.js Features (Built-in)
* Mention to the user that the blog features an interactive 3D plugin (`pg_bio_sync.js`).
* The user can double-click either 3D viewer to lock their cameras together (synchronized rotation/tilt).
* The user can click any fragment on one protein, and it will automatically highlight the matching residue (across all complex chains!) on the opposite protein by swapping colors.
