#!/usr/bin/env python3

# Usage: 
#     generate-authors.py > AUTHORS


# Authors
authors = [
    "Willem Deconinck",
    "Pedro Maciel",
    "Tiago Quintino"
]

def real_name(name):
    # Because some contributors use different names or aliases
    alias = {
        "cosunae"               : "Carlos Osuna",
        "carlos osuna"          : "Carlos Osuna",
        "mengaldo"              : "Gianmarco Mengaldo",
        "svahl991"              : "Steven Vahl",
        "benjaminmenetrier"     : "Benjamin Menetrier",
        "danholdaway"           : "Daniel Holdaway",
        "MarekWlasak"           : "Marek Wlasak",
        "odlomax"               : "Oliver Lomax",
        "MO-marcomilan"         : "Marco Milan",
        "twsearle"              : "Toby Searle",
        "Dusan Figala"          : "Dušan Figala",
    }
    if name in alias:
        return alias[name]
    else:
        return name

def contributors():
    # Contributors sorted by number of commits
    from subprocess import check_output
    command=['git', 'shortlog', '-n', '-s']
    output=check_output(command, universal_newlines=True).strip()
    contributors_dict = {}
    for line in output.split('\n'):
        [commits,name] = line.strip().split('\t',1)
        name = real_name(name)
        if name in authors:
            continue
        commits = int(commits)
        if name in contributors_dict:
            commits += contributors_dict[name]
        contributors_dict[name] = commits
    return [item[0] for item in sorted(contributors_dict.items(), key = lambda x: x[1], reverse = True)]

print("Authors")
print("=======\n")
for name in authors:
  print("- "+name)

print("")
print("Thanks for contributions from")
print("=============================\n")

for name in contributors():
  print("- "+name)
