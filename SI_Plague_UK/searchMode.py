

import os

os.system('plom bootstrap theta.json design_lhs.json --run --force -w "00:00:10"')
os.system('plom reduce design_lhs.json -o mle.json')
os.system('plom pipe mle.json | ./ksimplex --prior -M 1000')
for i in range(10):
    os.system('plom pipe mle.json -T | ./ksimplex --prior -M 1000')
os.system('plom bootstrap mle.json design_sample.json --run --force')
os.system('plom reduce design_sample.json -o mle.json')
