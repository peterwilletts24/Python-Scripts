def LinePlotEMBRACEExperimentID(experiment_id):
    '''
    Takes EMBRACE experiment_id as input and returns line colour, thickness and style
    '''

    NUM_COLOURS = 17
    from matplotlib import cm
    cmap=cm.get_cmap(cm.Set1, NUM_COLOURS)

    if (experiment_id=='djznw'):
        print experiment_id
        colour = 'black'
        linewidth=1.2
        linestylez='--'
    if (experiment_id=='djzny'):
        print experiment_id
        colour = cmap(1.*3/NUM_COLOURS)
        linewidth=0.5
        linestylez='--'
    if ((experiment_id=='djznq') or (experiment_id=='dkjxq')):
        print experiment_id
        colour = cmap(1.*11/NUM_COLOURS)
        linewidth=0.8
        if (experiment_id=='djznq'):
            linestylez='--'
        if (experiment_id=='dkjxq'):
            linestylez=':'
              
    if ((experiment_id=='dklzq') or (experiment_id=='dklwu')):
            print experiment_id
            colour = cmap(1.*7/NUM_COLOURS)
            linewidth=1
            if (experiment_id=='dklzq'):
                linestylez='--'
            if (experiment_id=='dklwu'):
                linestylez='-'
    if ((experiment_id=='dklyu') or (experiment_id=='dkmbq')):
            print experiment_id
            colour = cmap(1.*9/NUM_COLOURS)
            linewidth=1.3
            if (experiment_id=='dkmbq'):
                linestylez='--'
            if (experiment_id=='dklyu'):
                linestylez='-'
    if (experiment_id=='djzns'):
                print experiment_id
                colour = cmap(1.*11/NUM_COLOURS)
                linewidth=1.6
                linestylez='-'
    if ((experiment_id=='dkbhu') or (experiment_id=='dkhgu')):
             print experiment_id
             colour = cmap(1.*13/NUM_COLOURS)
             linewidth=1.9
             if (experiment_id=='dkbhu'):
                 linestylez='-'
             if (experiment_id=='dkhgu'):
                 linestylez=':'
    if (experiment_id=='djznu'):
             print experiment_id
             colour = cmap(1.*15/NUM_COLOURS)
             linewidth=2.
             linestylez='-'        

    if (experiment_id=='dkmgw'):
             print experiment_id
             colour = cmap(1.*16/NUM_COLOURS)
             linewidth=1.3
             linestylez='-'

    return colour, linewidth, linestylez
