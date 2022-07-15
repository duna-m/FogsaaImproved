from oldMain import *
import matplotlib.pyplot as plt
import numpy as np

if __name__ == '__main__':
    args = sys.argv[1:]
    '''args[0] shouolld be evyat0.txt'''

    probabilities = np.array([[0,0,0],[0.01,0.01,0.01],[0.03,0.03,0.03],[0.06,0.06,0.06],[0.09,0.09,0.09],[0.12,0.12,0.12],[0.15,0.15,0.15],[0.18,0.18,0.18],[0.21,0.21,0.21],[0.24,0.24,0.24],[0.27,0.27,0.27],[0.3,0.3,0.3],[0.33,0.33,0.33]])
    wunsch_results = np.empty((0, 4), float)
    fogsaa_results = np.empty((0, 4), float)



    for i in range(len(probabilities)):
        evyat = "evyat" + str(i) + ".txt"
        similarity = 1 - (probabilities[i][0] + probabilities[i][1] + probabilities[i][2])
        wunsch_avg_time, wunsch_scs, wunsch_lcs, wunsch_lev, fogsaa_avg_time, fogsaa_scs, fogsaa_lcs, fogsaa_lev = \
            old_main([evyat, -(probabilities[i][0]*100), -(probabilities[i][1]*100), -(probabilities[i][2]*100), similarity])
        wunsch_results = np.append(wunsch_results, [[wunsch_avg_time, wunsch_scs, wunsch_lcs, wunsch_lev]], axis=0)
        fogsaa_results = np.append(fogsaa_results, [[fogsaa_avg_time, fogsaa_scs, fogsaa_lcs, fogsaa_lev]], axis=0)

    #wunsch_results = np.append(wunsch_results, np.array([[0,0,0,0],[1,1,1,1],[2,2,2,2],[3,3,3,3],[4,4,4,4],[5,5,5,5]]),axis = 0)
    #fogsaa_results = np.append(fogsaa_results, np.array([[0,0,0,0],[3,3,3,3],[6,6,6,6],[9,9,9,9],[12,12,12,12],[15,15,15,15]]),axis = 0)

    fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4)

    l1 = ax1.plot(probabilities[:, 0], wunsch_results[:, 0],marker='o', label='Needleman-Wunsch', color="red")
    l2 = ax1.plot(probabilities[:, 0], fogsaa_results[:, 0],marker='x', label='Fogsaa', color="blue")
    ax1.title.set_text('Avg Time')
    #ax1.legend(loc='upper left')
    ax1.set(xlabel='insertion probability')


    l3 = ax2.plot(probabilities[:, 0], wunsch_results[:, 1],marker='o', label='Needleman-Wunsch', color="red")
    l4 = ax2.plot(probabilities[:, 0], fogsaa_results[:, 1],marker='x', label='Fogsaa', color="blue")
    ax2.title.set_text('SCS metric')
    #ax2.legend(loc='upper left')
    ax2.set(xlabel='insertion probability')




    l5 = ax3.plot(probabilities[:, 0], wunsch_results[:, 2],marker='o', label='Needleman-Wunsch', color="red")
    l6 = ax3.plot(probabilities[:, 0], fogsaa_results[:, 2],marker='x', label='Fogsaa', color="blue")
    ax3.title.set_text('LCS metric')
    #ax3.legend(loc='upper left')
    ax3.set(xlabel='insertion probability')




    l7 = ax4.plot(probabilities[:, 0], wunsch_results[:, 3],marker='o', label='Needleman-Wunsch', color="red")
    l8 = ax4.plot(probabilities[:, 0], fogsaa_results[:, 3],marker='x', label='Fogsaa', color="blue")
    ax4.title.set_text('Levenstein metric')
    #ax4.legend(loc='upper left')
    ax4.set(xlabel='insertion probability')



    fig.suptitle('metrics as a function of insertion probability')
    fig.legend([l1, l2],  # The line objects
               labels=["needleman-wunsch", "fogsaa"],  # The labels for each line
               loc="center right",  # Position of legend
               borderaxespad=0.1,  # Small spacing around legend box
               title="Legend Title"  # Title for the legend
               )
    plt.show()









