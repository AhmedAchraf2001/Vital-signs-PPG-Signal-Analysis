def qppgfast_beat_detector(sig, fs):
    onsets = qppg_fast(sig,fs)
    onsets = tidy_beats(onsets)
    peaks = pulse_peaks_from_onsets(sig, onsets)
    return (peaks, onsets)


def qppg_fast(data, fs, form = 0 , to=7500):   

    # PC: if no start and end point are specified, then analyse the entire signal.
    import numpy as np
    import math
    fs = 125
    form = 0
    to = len(data)
    # PC: Setup
    ppgOnsets = []
    beat_n = 0

    sps = fs      # Sampling Frequency

    BUFLN = 4096            # /* must be a power of 2, see slpsamp() */
    EYE_CLS = 0.34          # /* eye-closing period is set to 0.34 sec (340 ms) for PPG */ 
    LPERIOD  = sps*8        # /* learning period is the first LPERIOD samples */
    SLPW = 0.17             # /* Slope width (170ms) for PPG */                        
    NDP = 2.5               # /* adjust threshold if no pulse found in NDP seconds */
    TmDEF = 5               # /* minimum threshold value (default) */
    Tm = TmDEF 

    INVALID_DATA = -32768
    c=0

    if data[0] <= INVALID_DATA+10 :
        data[0] = np.mean(data)
   
    inv = np.where( data <= INVALID_DATA + 10 )
    for i in range(1, len(inv)):
        data[inv[i]] = data[inv[i] - 1]

    # re-scale data to ~ +/- 2000
    if len(data) < 5*60*sps :
        data = np.divide(( data - min(data) ), ( max(data) - min(data) ))* 4000 - 2000
    else :
        n=1 
        max_data = []
        min_data = []
        for i in range(1, len(data), 5*60*sps):
            max_data.append( max( data[ i:min( i+5*60*sps-1,len(data) ) ] ) )
            min_data.append( min( data[ i:min( i+5*60*sps-1,len(data) ) ] ) ) 
            n=n+1 

        data = np.divide( ( data - np.median(min_data) ), ( np.median(max_data) - np.median(min_data) ))* 4000-2000 
    

    EyeClosing = math.ceil(sps * EYE_CLS)    # /* set eye-closing period */
    ExpectPeriod = math.ceil(sps * NDP)      # /* maximum expected RR interval */
    SLPwindow = round(sps * SLPW)        # /* slope window size */
    timer = 0 

    # PC: Setup variables
    ebuf = np.zeros((1, BUFLN))
    ebuf = ebuf[0]
    lbuf = np.zeros((1, BUFLN))
    lbuf = lbuf[0] 

    if form > BUFLN :
        tt_2 = form - BUFLN 
    else :
        tt_2 = 0 

    aet = 0 

    t1= 8*sps 
    t1 = t1 + form 
    T0 = 0 
    n = 0  
    
    # loop through the first 8 seconds of the signal
    for t in range(form, t1):
        temp,ebuf,lbuf,tt_2, aet = slpsamp(t,data,BUFLN,ebuf,lbuf,tt_2, aet,SLPwindow) 
        if temp > INVALID_DATA+10 :
            T0 = T0+temp 
            n=n+1 
            
    T0 = T0/n         # T0=T0/(t1-form) 
    Ta = 3 * T0

    learning = 1

    
    
    #    /* Main loop */

    t = form
    while t <= to :
        if (learning) :
            if (t > (form + LPERIOD)):
                learning = 0
                T1 = T0
                t = form        # /* start over */
            else :
                T1 = 2*T0

        temp,ebuf,lbuf,tt_2, aet = slpsamp(t,data,BUFLN,ebuf,lbuf,tt_2, aet,SLPwindow)

        if (temp > T1):    # /* found a possible ABP pulse near t */ 
            timer = 0                # /* used for counting the time after previous ABP pulse */
            maxd = temp
            mind = maxd
            tmax = t
            for tt in range(t + 1, t + EyeClosing):
                [temp2 ,ebuf,lbuf,tt_2, aet] = slpsamp(tt,data,BUFLN,ebuf,lbuf,tt_2, aet,SLPwindow)
                if temp2 > maxd :
                    maxd=temp2
                    tmax=tt
            if maxd == temp:
                t=t+1
                continue
            
            for tt in range(tmax, int(t - EyeClosing / 2 ), -1) :
                [temp2 ,ebuf,lbuf,tt_2, aet] = slpsamp(tt,data,BUFLN,ebuf,lbuf,tt_2, aet,SLPwindow)
                if temp2 < mind:
                    mind = temp2
            if maxd > mind + 10 :
                onset = (maxd - mind)/100+2
                tpq = t - round(0.04*fs)
                maxmin_2_3_threshold = (maxd - mind)*2.0/3
                for tt in range(tmax, int(t-EyeClosing/2), -1):
                    [temp2, ebuf,lbuf,tt_2, aet] = slpsamp(tt,data,BUFLN,ebuf,lbuf,tt_2, aet,SLPwindow)
                    if temp2 < maxmin_2_3_threshold:
                        break
                
                for tt in range(tt, int(t - EyeClosing / 2 + round(0.024*fs) -1 ), -1):
                    [temp2 ,ebuf,lbuf,tt_2, aet] = slpsamp(tt,data,BUFLN,ebuf,lbuf,tt_2, aet,SLPwindow)
                    [temp3 ,ebuf,lbuf,tt_2, aet] = slpsamp(tt-round(0.024*fs),data,BUFLN,ebuf,lbuf,tt_2, aet,SLPwindow)
                    if temp2-temp3 < onset :
                        tpq = tt- round(0.016*fs)
                        break
                
                # find valley form the original data around 0.25s of tpq 
                valley_v = round(tpq)
                for valley_i in range(round(max(2,tpq-round(0.20*fs))), round(min(tpq+round(0.05*fs),len(data)-1)+1)):
                    # If vally is too low, it cannot serve as an index, so move to the next time.
                    if valley_v <= 0 :
                        t += 1
                        continue
                
                    if data[valley_v] > data[valley_i] and data[valley_i] <= data[valley_i-1] and data[valley_i] <= data[valley_i+1]:
                        valley_v = valley_i

                
                if not(learning) :
                    # If we are looking for the first peak
                    if beat_n == 0 :
                        # If the proposed peak index > 0
                        if round(valley_v) > 0 :
                            ppgOnsets.insert(beat_n, round(valley_v))
                            beat_n += 1
                    else :
                        # Check if rounded valley_v is greater than the prior beat index
                        if round(valley_v) > ppgOnsets[beat_n-1] :
                            ppgOnsets.insert(beat_n, round(valley_v))
                            beat_n += 1
                        

                # /* Adjust thresholds */
                Ta = Ta + (maxd - Ta)/10
                T1 = Ta / 3

                # /* Lock out further detections during the eye-closing period */
                t = tpq + EyeClosing
            
        else :
            c+=1
            if not(learning) :
                #  Once past the learning period, decrease threshold if no pulse was detected recently. 
                timer = timer + 1
                if timer > ExpectPeriod and Ta > Tm : 
                    Ta = Ta - 1
                    T1 = Ta / 3

        t=t+1

    # Discard first beat because algorithm always finds first minimum value, so trace-back logic will find a fir
    return ppgOnsets


def slpsamp(t, data, BUFLN, ebuf, lbuf, tt_2, aet, SLPwindow):
    import numpy as np
    while (t > tt_2):
        prevVal = 0
        if (tt_2 > 0) and (tt_2 -1 > 0) and (tt_2 < len(data)) and (tt_2 -1 < len(data)):
            val2 = data[tt_2 - 2]
            val1 = data[tt_2 - 1]
        else :
            val2 = prevVal
            val1 = val2
        prevVal = val2
        dy =  val1 - val2

        if dy < 0 :
            dy = 0
        tt_2 = tt_2 + 1 
        M = np.round(np.remainder(tt_2, (BUFLN - 1)) + 1)
        et = dy
        ebuf[M-1] = et
        aet = 0
        for i in range(SLPwindow):
            p = M - i 
            if p < 0 :
                p += BUFLN
            aet += ebuf[p-1]
        lbuf[M-1] = aet
    M3 = np.round(np.remainder(t, (BUFLN - 1)) + 1)
    beat1 = lbuf[M3-1]
    return (beat1, ebuf, lbuf, tt_2, aet)