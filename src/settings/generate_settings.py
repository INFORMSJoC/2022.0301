import json

with open("default_settings.json", "r") as read_file:
    data = json.load(read_file)

for scoring in ["MinFunction", "MaxFunction", "ProductScoring", "WeightedSum", "WeightedProduct"]:
    for st_branching in [1, 4, 8, 10, 12, 16]:
        for br_point in [0.2, 0.5, 0.7]:
            data["scoring_parameter"] = scoring
            data["strong_branching"] = st_branching
            data["branching_point"] = br_point
            with open("{}-SB_{}-BP_{}.json".format(scoring, st_branching, br_point), "w") as write_file:
                json.dump(data, write_file)
