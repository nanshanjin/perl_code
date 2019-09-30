import os

for i in os.listdir("/home/Account/yangns/project/metagenomics/JZ201907311146_Bacteria_Fungi/02-kraken2"):
    if ".report.xls" in i:
        print i
