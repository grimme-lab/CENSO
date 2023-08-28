import json

def main():
    with open("./censo_dfa_settings.json", "r") as f:
        data = json.load(f)

    newdata = data.copy()
    newdata.pop("relay_functionals")

    for k, v in newdata.get("functionals").items():
        if "wb97" not in k and "dsd-blyp" not in k and "3c" not in k:
            func = k.split()
            for prog in ["tm", "orca"]:
                if v.get(prog) == "":
                    v.pop(prog)
            
            if func[0] not in newdata.get("functionals").keys():
                newdata.get("functionals")[func[0]] = {"type": v.get("type"), "disp": []}
                for prog in ["tm", "orca"]:
                    if v.get(prog):
                        newdata.get("functionals")[func[0]][prog] = v.get(prog)
            
            newdata.get("functionals")[func[0]]["disp"].append(func[1])
        elif "3c" in k:
            newdata.get("functionals")[k] = v


    with open("./censo_dfa_settings.json", "w") as f:
        json.dump(data, f, indent=4)


if __name__ == "__main__":
    main()
