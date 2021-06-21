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
            "part4": None,
            "part5": None,
        },
        previous_part_info={
            "part0": 0.0,
            "part1_firstsort": 0.0,
            "part1": 0.0,
            "part2_opt": 0.0,
            "part2": 0.0,
            "part3": 0.0,
            "part4": 0.0,
            "part5": 0.0,
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
            "part4": None,
            "part5": None,
        },
        supporting_info={
            "part0": {
                "Energy": None,
                "Energy_settings": None,
                "G_mRRHO": None,
                "G_solv": None,
                "Geometry": None,
                "Threshold": None,
                "main QM code": None,
            },
            "part1": {
                "Energy": None,
                "Energy_settings": None,
                "G_mRRHO": None,
                "G_solv": None,
                "Geometry": None,
                "Threshold": None,
                "main QM code": None,
            },
            "part2": {
                "Energy": None,
                "Energy_settings": None,
                "G_mRRHO": None,
                "G_solv": None,
                "Geometry": None,
                "Threshold": None,
                "main QM code": None,
            },
            "part3": {
                "Energy": None,
                "Energy_settings": None,
                "G_mRRHO": None,
                "G_solv": None,
                "Geometry": None,
                "Threshold": None,
                "main QM code": None,
            },
            "part4": {
                "Energy": None,
                "Energy_settings": None,
                "G_mRRHO": None,
                "G_solv": None,
                "Geometry": None,
                "Threshold": None,
                "main QM code": None,
            },
            "part5": {
                "Energy": None,
                "Energy_settings": None,
                "G_mRRHO": None,
                "G_solv": None,
                "Geometry": None,
                "Threshold": None,
                "main QM code": None,
            },
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
        self.previous_part_info = previous_part_info
        self.avGcorrection = avGcorrection
        self.comment = comment
        self.bestconf = bestconf
        self.nconfs_per_part = nconfs_per_part
        self.si = supporting_info
