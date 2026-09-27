import os

filepath = "../ideia.me/_posts/2026-09-27-decoding-ai-vector-distances.md"
with open(filepath, "r") as f:
    content = f.read()

new_content = """---

> **A Quick Jargon Check: What is a "Bait"?**
> In biology, the term "Bait" comes from a very literal fishing analogy. When you go fishing, you put something you already know and possess (the bait) onto a hook, cast it into a massive, dark ocean, and wait to see what unknown things "bite" or stick to it. In computational biology, our "Bait" is the 1024-dimensional vector of a well-documented protein. We cast it into the database to see which unknown proteins have the exact same shape and "stick" to it!

---

## Visualizing the 0.080 Distance 

Don't just trust the math—trust your eyes. Below is a 3D visualization comparing a known Thiamine-binding "Bait" protein against our top discovery (Distance: `0.080`). 

When the Cosine Distance drops this low, you can visually see that the physical architectures—the folding of the Alpha-helices and Beta-sheets—are functionally identical, despite coming from completely different branches of the Tree of Life.

<!-- Include the 3Dmol.js Library -->
<script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>

<div style="display: flex; gap: 20px; justify-content: center; flex-wrap: wrap; margin-top: 20px;">
    <!-- Known Bait Protein -->
    <div style="text-align: center;">
        <h4>Known Bait (Thiamine-binding)</h4>
        <div style="height: 400px; width: 350px; position: relative; border: 1px solid #ccc; border-radius: 8px;" 
             class="viewer_3Dmoljs" 
             data-href="/assets/pdb/thiamine_bait.pdb" 
             data-backgroundcolor="0x1e1e1e" 
             data-style="cartoon:color=orange">
        </div>
    </div>

    <!-- Orphan Discovery -->
    <div style="text-align: center;">
        <h4>Our Discovery (Distance: 0.080)</h4>
        <div style="height: 400px; width: 350px; position: relative; border: 1px solid #ccc; border-radius: 8px;" 
             class="viewer_3Dmoljs" 
             data-href="/assets/pdb/thiamine_orphan.pdb" 
             data-backgroundcolor="0x1e1e1e" 
             data-style="cartoon:color=purple">
        </div>
    </div>
</div>

<p style="text-align: center; font-style: italic; margin-top: 10px;">
(Notice the deep, central pocket shared by both proteins where the Thiamine molecule physically docks!)
</p>

---

## Conclusion: Math is the New Microscope"""

old_anchor = "---\n\n## Conclusion: Math is the New Microscope"

if old_anchor in content:
    content = content.replace(old_anchor, new_content)
    with open(filepath, "w") as f:
        f.write(content)
    print("Successfully added Bait explanation and 3D preview.")
else:
    print("Could not find anchor to replace.")
