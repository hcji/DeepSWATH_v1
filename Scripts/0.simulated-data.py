# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 09:51:17 2020

@author: hcji
"""

import math
import random
import numpy as np
import matplotlib.pyplot as plt

def generate_single_profile(decoy = False):
    x = np.arange(100)
    if decoy:
        mu = [random.uniform(20, 40), random.uniform(60, 80)]
        mu = random.choice(mu)
    else:        
        mu = 50
    sigma = random.uniform(1,10)
    y = np.exp(-1*((x-mu)**2)/(2*(sigma**2)))/(math.sqrt(2*np.pi) * sigma)
    y /= max(y)
    # plt.plot(x, y)
    return y


def generate_overlap_profile(decoy = False):
    x = np.arange(100)
    y1 = generate_single_profile(decoy)
    for i in range(random.randint(1,3)):
        mu = [random.uniform(20, 40), random.uniform(60, 80)]
        mu = random.choice(mu)
        sigma = random.uniform(1, 10)
        y2 = np.exp(-1*((x-mu)**2)/(2*(sigma**2)))/(math.sqrt(2*np.pi) * sigma)
        y2 = y2 / max(y2) * random.uniform(0.3, 3)
        y1 += y2
    y1 = y1 / max(y1)
    # plt.plot(x, y)
    return y1
    
    
def add_noise(y):
    y1 = y
    y1 = y1 * np.array([ random.uniform(0.8, 1.2) for i in range(len(y))])
    y1 = y1 + np.array([ random.uniform(-0.3, 0.3) for i in range(len(y))])
    y1 /= max(y1)
    # plt.plot(y1)
    return y1


if __name__ == '__main__':
    
    from tqdm import tqdm
    
    simu_precursor_eics = []
    simu_fragment_eics = []
    simu_decoy_eics = []    
    
    for i in tqdm(range(72214)):
        
        if random.choice([True, False]):
            x = generate_single_profile()
        else:
            x = generate_overlap_profile()
            
        if random.choice([True, False]):        
            y1 = generate_single_profile()
        else:
            y1 = generate_overlap_profile()
            
        if random.choice([True, False]):
            y2 = generate_single_profile(decoy = True)
        else:
            y2 = generate_overlap_profile(decoy = True)
 
        x = add_noise(x)
        y1 = add_noise(y1)
        y2 = add_noise(y2)
        '''
        plt.plot(x)
        plt.plot(y1)
        plt.plot(y2)
        '''
        simu_precursor_eics.append(x)
        simu_fragment_eics.append(y1)
        simu_decoy_eics.append(y2)
    
    simu_precursor_eics = np.asarray(simu_precursor_eics)
    simu_fragment_eics = np.asarray(simu_fragment_eics)
    simu_decoy_eics = np.asarray(simu_decoy_eics)
    np.save('Data/simu_precursor_eics.npy', simu_precursor_eics)
    np.save('Data/simu_fragment_eics.npy', simu_fragment_eics)
    np.save('Data/simu_decoy_eics.npy', simu_decoy_eics)
    
    