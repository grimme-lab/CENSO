from typing import Union
import yaml
from dataclasses import dataclass

@dataclass
class Setting:
    type: type
    part: str
    name: str
    value: Union[int, float, str, list, bool]
    
    def to_yaml(self) -> dict:
        return {
            "type": self.type.__name__, 
            "part": self.part, 
            "name": self.name, 
            "value": self.value
        }
        
def main():
    settings = []
    
    settings.append(Setting(int, "prescreening", "nconf", 4))
    settings.append(Setting(list, "optrot", "freq_or", [12.3, 1.5]))
    
    with open("test", "w") as file:
        for setting in settings:
            yaml.safe_dump(setting.to_yaml(), file)
            
        
if __name__ == "__main__":
    main()