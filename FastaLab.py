#!/usr/bin/env python

#reads in fastq files with corresponding quality scores in a list
with open('Fastafile.fq', 'r') as f:
    lines = f.readlines()

#basically loops through lines in the file
for i in range(0, len(lines), 4):
    seq_id = lines[i].rstrip()  #extracts sequence ID
    seq = lines[i+1].rstrip()     #extracts sequence
    scores = lines[i+3].rstrip()    #extracts quality scores

    window_size = 10     #window size
    threshold = 15    #quality score less than 15
    low_quality_regions = []    #putting scores in a array
    new_seq = ''     #processing sequence
    new_scores = ''   #processing scores

    for j in range(0, len(seq)-window_size+1, window_size):
        window_scores = [ord(c)-33 for c in scores[j:j+window_size]]  #calculating the average quality score for the current window
        avg_score = sum(window_scores) / len(window_scores)

        if avg_score < threshold:
            low_quality_regions.append((j, j+window_size)) #if the avg score is below the threshold, marks as low quality
        else:
            new_seq += seq[j:j+window_size]   #if the avg quality score is above the threshold, appends the current window to the new sequence
            new_scores += scores[j:j+window_size]

    new_seq_final = ''  #removes stretches of DNA with low quality scores
    new_scores_final = ''
    i = 0
    while i < len(new_seq):
        is_low_quality = False  #this checks if the current position is part of a low quality region
        for start, end in low_quality_regions:
            if i >= start and i < end:
                is_low_quality = True
                i = end
                break

        if not is_low_quality:
            new_seq_final += new_seq[i]  # if its not part of low quality regions, then append to final sequence
            new_scores_final += new_scores[i]
            i += 1

    print(seq_id)  #prints outputs
    print(new_seq_final)
    print('+')
    print(new_scores_final)