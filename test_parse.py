import csv

file = csv.reader(
       filter(lambda line: line[0]!='\n',
              open('test_data.csv')),
       delimiter=" ")
for row in file:
    values = list(map(float, row))
    print(values)
