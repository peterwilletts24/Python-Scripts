'Converts model run names e.g djznu to sub title for plots'


def main(experiment_id): 
    djznu = ('1.5km', '4500 km', '10s', 'L118, 78km lid', 'Explicit 3D SMAG', 'djznu', '1.5kmExp')
    dkbhu = ('2.2km', '4500 km', '10s', 'L118, 78km lid', 'Explicit 3D SMAG', 'dkbhu', '2.2kmExp')
    djzns = ('4km', '4500 km', '10s', 'L118, 78km lid', 'Explicit 3D SMAG', 'djzns', '4kmExp')
    djznq = ('24km',	'4500 km',	'600s',	'Global NWP set L70, 80 km lid', 	'1DBL + conv param', 'djznq', '24kmExp')
    djzny = ('120km',	'4500 km',	'1200s',	'Global NWP set L70, 80 km lid', 	'1DBL + conv param', 'djzny', '120kmExp')
    djznw = ('Driving Global', '4500 km', '1200s?', 'Global NWP set L70, 80 km lid', '1DBL + conv param', 'djznw', 'DrivingLAM')
						
    dkhgu = ('2.2km Big', 'Big', '', '', '', 'dkhgu', '2.2kmBig')			
    dkjxq = ('24km Big','Big', '', '', '', 'dkjxq', '24kmBig')
						
    dklyu = ('8km','4500 km' ,'10s', 'L118, 78km lid', 'Explicit 3D SMAG', 'dklyu', '8kmExp')
    dkmbq = ('8km', '4500 km', '300s?', '', '1DBL + conv param', 'dkmbq', '8kmPar')
						
    dklwu = ('12km', '', '10s', 'L118, 78km lid', 'Explicit 3D SMAG', 'dklwu', '12kmExp')
    dklzq = ('12km', '', '300s?', '1DBL + conv param', 'dklzq', '12kmPar')

    experiment_ids = [djznu, dkbhu, djzns, djznq, djzny, djznw, dkhgu, dkjxq, dklyu, dkmbq, dklwu, dklzq ]


    for ex in experiment_ids:

        if (experiment_id==ex[5]):
            #mod_sub_title = ex


             return ex

# 'Run: %s, Resolution: %s, Domain: %s, Timestep: %s, Vert. Levels %s, Conv. scheme %s' % (title[-1], title[0],  title[1], title[2], title[3], title[4])

   

if __name__ == '__main__':
    main()
