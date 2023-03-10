"""
Contains info regarding the ensemble such as energies etc., best conf per part,
number of confs per part, timings (remove)
metadata class
"""

# TODO - should contain all information calculated results for all conformers
# identify data by conformer id
# fairly big/highly functional class
class EnsembleData:
    def __init__(
        self,
        id="ensemble_info",
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
        # TODO - reduce copy/paste code
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
        self.avGcorrection = avGcorrection
        self.comment = comment
        self.bestconf = bestconf
        self.nconfs_per_part = nconfs_per_part
        self.si = supporting_info
