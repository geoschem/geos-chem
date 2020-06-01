#!/usr/bin/env python

import yaml

file_name = "species_database.yml"
species_db = yaml.load(open(file_name), Loader=yaml.FullLoader)
print("species db")
print(type(species_db))

print(species_db.keys())

file_name = "species_database_tomas.yml"
species_db = yaml.load(open(file_name), Loader=yaml.FullLoader)
print("species db TOMAS")
print(species_db.keys())

file_name = "species_database_apm.yml"
species_db = yaml.load(open(file_name), Loader=yaml.FullLoader)
print("species db APM")
print(species_db.keys())
