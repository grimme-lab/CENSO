class EnsembleData:
    def __init__(
        self,
        id="ensemble_info",
        filename=None,
        part_info={
            "part0": None,
            "part1_firstsort": None,
            "part1": None,
            "part2_opt": None,
            "part2": None,
            "part3": None,
        },
        avGcorrection=None,
        comment=None,
        bestconf={"part0": None, "part1": None, "part2": None, "part3": None},
        nconfs_per_part={
            "starting": None,
            "part0": None,
            "part1_firstsort": None,
            "part1": None,
            "part2_opt": None,
            "part2": None,
            "part3": None,
        },
    ):
        """
       ensemble_data: Creates an object where data
                       concerning the entire ensemble is stored.
        Input:
        filename = e.g. crest_conformers.xyz
        part_info --> time passed to calculate part
        avGcorrection --> information of higher lying conformers
        bestconf --> id of best conf per part
        nconfs_per_part --> how many confs have been evaluated in each part

        """
        if avGcorrection is None:
            avGcorrection = {}
        if comment is None:
            comment = []
        self.id = id
        self.filename = filename
        self.part_info = part_info
        self.avGcorrection = avGcorrection
        self.comment = comment
        self.bestconf = bestconf
        self.nconfs_per_part = nconfs_per_part
