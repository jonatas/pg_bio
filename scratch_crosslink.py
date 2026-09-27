import os

posts_dir = "../ideia.me/_posts/"
files = [
    "2026-09-25-why-proteins-are-just-high-dimensional-vectors.md",
    "2026-09-25-writing-a-hybrid-database-operator-for-biology.md",
    "2026-09-25-indexing-3d-space-like-a-video-game.md",
    "2026-09-25-mining-the-dark-proteome-legacy-code-with-missing-docs.md"
]

slugs = [
    "/why-proteins-are-just-high-dimensional-vectors",
    "/writing-a-hybrid-database-operator-for-biology",
    "/indexing-3d-space-like-a-video-game",
    "/mining-the-dark-proteome-legacy-code-with-missing-docs"
]

titles = [
    "Part 1: Why Proteins are just High-Dimensional Vectors",
    "Part 2: Writing a Hybrid Database Operator for Biology (<~>)",
    "Part 3: Indexing 3D Space like a Video Game",
    "Part 4: Mining the Dark Proteome: Legacy Code with Missing Docs"
]

def generate_header(index):
    lines = [
        "> **The pg_bio Series**",
        "> This post is part of a 4-part series on building a bioinformatics engine natively in PostgreSQL.",
        "> "
    ]
    for i in range(4):
        if i == index:
            lines.append(f"> * **{titles[i]}** (You are here)")
        else:
            lines.append(f"> * [{titles[i]}]({slugs[i]})")
    lines.append("\n---\n")
    return "\n".join(lines)

def generate_footer(index):
    if index < 3:
        return f"\n---\n\n**Next up in the series:** [{titles[index+1]}]({slugs[index+1]})"
    else:
        return f"\n---\n\n*This concludes the pg_bio series! Check out [Part 1]({slugs[0]}) if you missed how it all started.*"

for i, f in enumerate(files):
    filepath = os.path.join(posts_dir, f)
    with open(filepath, "r") as file:
        content = file.read()
    
    # Split frontmatter
    parts = content.split("---", 2)
    if len(parts) == 3:
        frontmatter = parts[1]
        body = parts[2]
        
        # Add a preview separator (<!--more-->) after the first couple of paragraphs
        # Find the second double newline
        preview_split = body.split("\n\n", 3)
        if len(preview_split) >= 4:
            body = f"{preview_split[0]}\n\n{preview_split[1]}\n\n<!--more-->\n\n{preview_split[2]}\n\n{preview_split[3]}"

        new_body = "\n" + generate_header(i) + body + generate_footer(i)
        
        # Write back
        new_content = f"---{frontmatter}---{new_body}"
        with open(filepath, "w") as file:
            file.write(new_content)
        print(f"Updated {f}")
