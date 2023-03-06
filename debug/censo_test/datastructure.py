

class MoleculeData:
    """
    MoleculeData contains an identifier of the conformer (id) and geometry information
    in order to keep the object small, since it has to be pickled for multiprocessing.Queue
    """

    def __init__(self, id: int, xyz):
        """"""