import os

# Fix config
config_path = "../ideia.me/_config.yml"
with open(config_path, "a") as f:
    f.write("\ncategories_path: /categories.html\n")

# Fix cross links in the posts
posts_dir = "../ideia.me/_posts/"
files = [
    "2026-09-25-why-proteins-are-just-high-dimensional-vectors.md",
    "2026-09-25-writing-a-hybrid-database-operator-for-biology.md",
    "2026-09-25-indexing-3d-space-like-a-video-game.md",
    "2026-09-25-mining-the-dark-proteome-legacy-code-with-missing-docs.md"
]

for file in files:
    filepath = os.path.join(posts_dir, file)
    with open(filepath, "r") as f:
        content = f.read()
    
    # Replace the extensionless paths with .html
    content = content.replace("(/why-proteins-are-just-high-dimensional-vectors)", "(/why-proteins-are-just-high-dimensional-vectors.html)")
    content = content.replace("(/writing-a-hybrid-database-operator-for-biology)", "(/writing-a-hybrid-database-operator-for-biology.html)")
    content = content.replace("(/indexing-3d-space-like-a-video-game)", "(/indexing-3d-space-like-a-video-game.html)")
    content = content.replace("(/mining-the-dark-proteome-legacy-code-with-missing-docs)", "(/mining-the-dark-proteome-legacy-code-with-missing-docs.html)")

    with open(filepath, "w") as f:
        f.write(content)

