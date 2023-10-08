import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def rolling_average(num_avg, time, mag): 
    num_data = len(time) - num_avg
    avg_mag = np.ndarray((num_data))
    avg_time = np.ndarray((num_data))
    for i in range(0,num_data):
        avg_mag[i] = np.average(mag[i:i+num_avg])
        avg_time[i] = np.average(time[i:i+num_avg])
    return(avg_mag, avg_time)

def sort_data(data):
    num_data = len(data)

    time = np.ndarray((num_data))
    mag = np.ndarray((num_data))
    err = np.ndarray((num_data))

    for i in range(0,num_data):
        time[i] = data[i][0]
        mag[i] = data[i][1]
        err[i] = data[i][2]
    return(time, mag, err)

def clean_data(err_lim, data):
    num_data=len(data)
    good_data = []
    for i in range(0, num_data):
        if data[i][2]<err_lim:
            good_data.append(i)
        else:
            continue
    
    cleaned = np.ndarray((len(good_data),3))
    for i in range(0, len(good_data)):
        cleaned[i]=data[good_data[i]]
    return(cleaned)

def time_norm(time,data):
    zero = time[data]
    for i in range(0,len(time)):
        time[i]=(time[i]-zero)*24
    return(time)

def mag_2_flux(mag, err):
    num_data = len(mag)
    flux = np.ndarray([num_data])
    flux_err = np.ndarray([num_data])
    for i in range(0,num_data):
        flux[i]=10**(mag[i]/-2.5)
        flux_err[i]=flux[i]*2.302585093*err[i]/-2.5
    return(flux, flux_err)

def norm_flux(flux, flux_err, norm):
    num_data = len(flux)
    norm_flux = np.ndarray([num_data])
    norm_flux_err = np.ndarray([num_data])
    for i in range(0,num_data):
        norm_flux[i]=flux[i]/norm
        norm_flux_err[i]=flux_err[i]/norm
    return(norm_flux, norm_flux_err)

def block_average(mag, err, time, block):
    half = int(block/2)
    end = int(len(time)-block/2)
    num = int(len(time)/block)
    block_mag = np.ndarray([num])
    block_time = np.ndarray([num])
    block_err = np.ndarray([num])
    for i in range(0,num):
        lower = int(i*block-block/2)
        upper = int(i*block+block/2)
        block_mag[i] = np.average(mag[lower:upper])
        block_time[i] = np.average(time[lower:upper])
        block_err[i] = np.sqrt(np.average(err[lower:upper]**2))
    return(block_time, block_mag, block_err)
    
        
        

    




    
