import os

filepath = "../ideia.me/resume.md"
with open(filepath, "r") as f:
    content = f.read()

old_intro = """I am a polymath and deeply curious builder. For over 20 years, I've worked across multiple disciplines—from hardcore database internals and AST parsing to generative art, financial algorithms, and now, structural biology. 

I believe the most exciting breakthroughs happen at the intersection of different fields. I love bridging domains, learning continuously, and sharing those discoveries through writing. Some examples of my explorations across different areas include:

* **Bioinformatics & Postgres:** [Indexing 3D Space like a Video Game](/indexing-3d-space-like-a-video-game.html)
* **AST Parsing & Tooling:** [Building a SQL Formatter with Fast](/building-a-sql-formatter-with-fast.html) 
* **Database Optimization:** [Benchmarking Ruby ORMs](/benchmarking-ruby-orms.html)
* **Creative Coding & Sonification:** [Mandala Playground - Interactive Sonification](/mandala-playground.html) """

new_intro = """I am a polymath, deeply curious builder, and dedicated [yogini](/yoga.html). For over 20 years, I've worked across multiple disciplines—from hardcore database internals and AST parsing to generative art, financial algorithms, and now, structural biology. 

I believe the most exciting breakthroughs happen at the intersection of different fields, connecting the technical with the physical and mindful. I love bridging domains, learning continuously, and sharing those discoveries through writing. Some examples of my explorations across different areas include:

* **Bioinformatics & Postgres:** [Indexing 3D Space like a Video Game](/indexing-3d-space-like-a-video-game.html)
* **AST Parsing & Tooling:** [Building a SQL Formatter with Fast](/building-a-sql-formatter-with-fast.html) 
* **Database Optimization:** [Benchmarking Ruby ORMs](/benchmarking-ruby-orms.html)
* **Creative Coding & Sonification:** [Mandala Playground - Interactive Sonification](/mandala-playground.html)
* **Mindfulness & Movement:** [My Yoga Practice](/yoga.html)"""

if old_intro in content:
    content = content.replace(old_intro, new_intro)
    with open(filepath, "w") as f:
        f.write(content)
    print("Successfully updated resume with yoga info.")
else:
    print("Could not find old string to replace.")
