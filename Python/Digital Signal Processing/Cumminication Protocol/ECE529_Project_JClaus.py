import math
import itertools as it
import numpy as np
import matplotlib.pyplot as plt
import random
import serial
import time
import csv

##Global Variables##
DC_shift = 3
Output_L = 2000
Output_Fs = 10000.0
Output_T = 1/Output_Fs
Input_L = 1000
Input_Fs = 10000.0
Input_T = 1/Input_Fs

##CSV Filenames##
key = 'key.csv'
device_data = 'device.csv'
d1_data = 'device1.csv'
d2_data = 'device2.csv'
d3_data = 'device3.csv'
d4_data = 'device4.csv'
d5_data = 'device5.csv'

##Raspberry Pi COM Port Settings##
#in_port = '/dev/ttyACM0'
#out_port = '/dev/ttyACM1'

##Laptop COM Port Settings##
in_port = 'COM2'
out_port = 'COM5'


#####Troubleshooting and Simulation Tools#####

## This function returns an array of random freq bits representative of
## the data array pointer's row number in binary with a
## size of n_freq = number of frequencies for the simulated
## device

def Create_TestFreqCall(N_Freq):
    bits = []
    for n in range(0,N_Freq):
        bits.append(random.randint(0, 1))
    return bits

def Simulate_SingleDevice():

    #Simulates a single device sending its data to a receiver
    bit_array = Create_TestFreqCall(50)
    freq_array = Build_FreqArray(bit_array, 100)
    out_signal = Build_Output(freq_array, Output_L, Output_T)
    comp_freqs = Find_Freq(out_signal, Input_Fs)

def Simulate_FiveDevices():

    #Simulates first device sending data 
    bit_array1 = Create_TestFreqCall(10)
    out_array1 = Build_FreqArray(bit_array1, 100)
    out_signal1 = Build_Output(out_array1, Output_L, Output_T)
    print('Device 1 O/P and data')
    print(out_array1)
    
    #Simulates second device receiving first device's data and
    #attaching its data to the signal 
    comp_freqs1 = Find_Freq(out_signal1, Input_Fs)
    print('I/P from Device 1')
    print(comp_freqs1)
    bit_array2 = Create_TestFreqCall(10)
    freq_array2 = Build_FreqArray(bit_array2, 1100)
    out_array2 = Append_FreqArray(comp_freqs1, freq_array2)
    out_signal2 = Build_Output(out_array2, Output_L, Output_T)
    print('Device 2 data')
    print(freq_array2)
    print('Device 2 O/P')
    print(out_array2)    

    #Simulates third device receiving second device's data and
    #attaching its data to the signal 
    comp_freqs2 = Find_Freq(out_signal2, Input_Fs)
    print('I/P from Device 2')
    print(comp_freqs2)
    bit_array3 = Create_TestFreqCall(10)
    freq_array3 = Build_FreqArray(bit_array3, 2100)
    out_array3 = Append_FreqArray(comp_freqs2, freq_array3)
    out_signal3 = Build_Output(out_array3, Output_L, Output_T)
    print('Device 3 data')
    print(freq_array3)
    print('Device 3 O/P')
    print(out_array3)

    #Simulates fourth device receiving third device's data and
    #attaching its data to the signal 
    comp_freqs3 = Find_Freq(out_signal3, Input_Fs)
    print('I/P from Device 3')
    print(comp_freqs3)
    bit_array4 = Create_TestFreqCall(10)
    freq_array4 = Build_FreqArray(bit_array4, 3100)
    out_array4 = Append_FreqArray(comp_freqs3, freq_array4)
    out_signal4 = Build_Output(out_array4, Output_L, Output_T)
    print('Device 4 data')
    print(freq_array4)
    print('Device 4 O/P')
    print(out_array4)

    #Simulates fifth device receiving fourth device's data and
    #attaching its data to the signal 
    comp_freqs4 = Find_Freq(out_signal4, Input_Fs)
    print('I/P from Device 4')
    print(comp_freqs4)
    bit_array5 = Create_TestFreqCall(10)
    freq_array5 = Build_FreqArray(bit_array5, 4100)
    out_array5 = Append_FreqArray(comp_freqs4, freq_array5)
    out_signal5 = Build_Output(out_array5, Output_L, Output_T)
    print('Device 5 data')
    print(freq_array5)
    print('Device 5 O/P')
    print(out_array5)

    #Simulates final monitoring device processing the data
    #from all five devices
    comp_freqs5 = Find_Freq(out_signal5, Input_Fs)
    print('I/P from Device 5')
    print(comp_freqs5)
    
#####TRANSMIT SIDE######

## This function converts the bit values represtative of the data row numbers (pointers)
## to a representative array of frequencies to be used by the Build_Output function


def Build_FreqArray(Bits, Min_Freq):

    # Defines the representative bits in the Bit Array as Freqencies in the Freq array
    Freq = []
    for i in range(0, len(Bits)):
        if (Bits[i] == 1.0):
            Freq.append((i*100.0) + Min_Freq)
            #print((i*100)+Min_Freq)
    #print(Freq)
    return Freq

## This function combines the component frequencies taken from an upstream device
## with the component frequencies for the current device to be used by the Build_Output
## function

def Append_FreqArray(Freq1, Freq2):

    for i in range(0, len(Freq2)):
        Freq1.append(Freq2[i])
    #print(Freq1)
    return Freq1
    


## This function builds an output array (OutputSignal) that is representative of the
## component frequencies from the frequency array (Freq_Array) and is scaled and biased
## so that the signal can be transmitted between 1 and 5 volts on 2 wires

def Build_Output(Out_FreqArray, L, T):

    # Initializes the Output Array
    OutputSignal = []
    for n in range(0, L-1):
        OutputSignal.append(0)

    # Builds the Output Array from the component frequencies and tracks the number of
    # frequencies (freq_count) for later scaling on output
    Freq_Count = 0
    for n in range(0, L-1):
        for i in range(0,len(Out_FreqArray)):
            OutputSignal[n] = OutputSignal[n] + math.sin(2*math.pi*Out_FreqArray[i]*n*T)
            Freq_Count = Freq_Count + 1
    plt.plot(OutputSignal)
    plt.show()

    # Scales and shifts the Output Array so that it can be transmitted as positive voltage
    # between 1 and 5 volts by tracking the additive properties of mutliple freqs and by
    # using a DC bias (global variable DC_shift)
    if (Freq_Count > 0):
        for n in range(0, L-1):
            OutputSignal[n] = (2 * L * OutputSignal[n] / Freq_Count) + DC_shift
    plt.plot(OutputSignal)
    plt.show()
    return OutputSignal

## This function sends the Output Signal array to be processed by the Arduino
## The Arduino continuously coverts the analog values in its buffer on a
## continuous loop to an 8 bit DAC input and stops its output to refresh its
## buffer with the values from the Output Array when receiving the "BeginBuffer"
## command. 

def Send_Signal(Signal):
    
    arduino_out = serial.Serial(out_port, 115200, timeout=.1)
    time.sleep(1) 
    arduino_out.write("BeginBuffer")
    for n in range(0,len(Signal)):
        arduino_out.write(Signal[n])
    arduino_out.write("EndBuffer")
    

#####RECEIVE SIDE#######

# This function reads the incoming signal from the Arduino and allocates the
# data to a list that is used for signal processing
def Read_Signal(L):

    arduino_in = serial.Serial(in_port, 115200, timeout=.1)
    data_count = 0
    signal_array = []
    while (data_count<=L):
        data_in = arduino_in.readline()[:-2] 
        if data_in:
            #print data
            signal_array.append(data_in)
            data_count = data_count + 1
    #print(signal_array)
    #plt.plot(signal_array)
    #plt.show()    
    return signal_array


#This funtion takes the Fourier Transform of the original signal and then computes the
#the single-sided spectrum from the double-sided spectrum P

def  Find_Freq(InputSignal,Fs):

    # Shifts the signal to be centered on the n axis
    for n in range(0, len(InputSignal)):
        InputSignal[n] = InputSignal[n] - DC_shift

    # Returns an fft of only the real single sided parts of the signal
    InputFreqs = np.fft.rfft(InputSignal)

    # Scales the spectrum to unity and tracks component freqs with the In_FreqArray
    for n in range(0, len(InputFreqs)):
        InputFreqs[n]=abs(InputFreqs[n]/len(InputFreqs))
    maxP = max(InputFreqs)
    plt.plot(InputFreqs)
    plt.show()    
    #print(maxP)
    In_FreqArray = []
    for n in range(0, len(InputFreqs)):
        InputFreqs[n]=abs(InputFreqs[n]/maxP)
        if (InputFreqs[n] > 0.6):
            In_FreqArray.append((1/2)*n*Fs/len(InputFreqs))
    plt.plot(InputFreqs)
    plt.show()
    #print(In_FreqArray)
    return In_FreqArray




#####Main Code#####



Simulate_SingleDevice()

#Simulate_FiveDevices()



