import json
import os
import pandas as pd

path = r"D:\python\3.MLST_SC\5.MLST\3.index_cal\3.candidate_genes_position_groups\2.position\candidate_cr01.json"
file = open(path)
f = json.load(file)
file.close()
# 在每个query中进行
query_list = []

result = {"name":[], "type":[], "start":[], "stop":[], "strand":[]}
querys = f["BlastOutput2"]
for query in querys:
    results = query["report"]["results"]["search"]
    # query_name = results["query_id"]
    query_name = results["query_title"]
    query_len = results["query_len"]
    hits = results["hits"]
    for hit in hits:
        hsp = hit["hsps"][0]
        result["name"].append(query_name)
        result["start"].append(min(hsp["hit_from"], hsp["hit_to"]))
        result["stop"].append(max(hsp["hit_from"], hsp["hit_to"]))
        if hsp["hit_strand"] == "Plus":
            result["strand"].append("+")
        else:
            result["strand"].append("-")
        result["type"].append("CDS")
        query_list.append(query_name)

df = pd.DataFrame(result)
print(df)
df.to_excel("position_groups_from_json.xlsx", index=False)

