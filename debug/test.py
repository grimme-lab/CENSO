from multiprocessing import Queue
        
from censo_test.settings import CensoSettings        

class Test:
    def __init__(self, q):
        self.q = q
        

    def put_in_q(self, item):
        self.q.put(item)
        
        
def main():
    q = Queue()
    
    settings = [CensoSettings()] * 100
    
    test = Test(q)
    for item in settings:
        test.put_in_q(item)
    test.put_in_q("STOP")

    ids = set()
    for i in iter(q.get, "STOP"):
        ids.add(id(i))
            
    print(ids)
        
if __name__ == "__main__":
    main()