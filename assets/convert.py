import json

def main():
    with open("./censo_solvents_db.json", "r") as f:
        dc = {}
        
        solvents = json.load(f)
        for solvent, definitions in solvents.items():
            dc[solvent] = definitions["DC"]

    with open("./solvents_dc.json", "w") as f:
        json.dump(dc, f, indent=4)


if __name__ == "__main__":
    main()
