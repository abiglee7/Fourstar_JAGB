import numpy as np
import pandas as pd
import os


def nir_trgb(galaxy, y_interc_J, mu_init,name, err_calcJ , y_intercep_H):

    """
    y_interc: expected TRGB magnitude
    mu_init: best guess for trgb distance modulus in T-space
    name: name of galaxy
    err_calc: the magnitude of the TRGB for which to estimate the average photometric uncertainty
    
    If there's an error about how it can't sigfig, then the muinit guess is probably off
    """
    
    for band in ["J", "H"]:
        
        if band=='J':
            x = np.arange(-1,1.6,.01)
            y_int = y_interc_J+.85; x_int = 0
            ones = np.ones_like(x)

            # rotate data around line of slope.85
            theta = np.arctan(.85)
            
            A = np.array([[np.cos(theta), np.sin(theta),0],[-np.sin(theta), np.cos(theta),0],[x_int-np.cos(theta)*x_int + y_int*np.sin(theta),\
                                y_int-x_int*np.sin(theta)-y_int*np.cos(theta),1]]).T
            
            ones_JK = np.ones(len(galaxy))
            
            
            JK_space = np.vstack((galaxy['J']-galaxy['K'], galaxy['J'], ones_JK))
            T_space = A.dot(JK_space)
            
            rotated_K = T_space[1] - T_space[0]
            
            # save to file for T-band TRGB plots if wanted
#             pd.DataFrame({'color':galaxy['J']-galaxy['K'], 'mag':T_space[1]}).to_csv('/Users/abigaillee/Documents/Fourstar project/Analysis/TRGB info/'+str(name)+'_tjtrgb.csv')
            
                   
            
            cat = pd.DataFrame({'I':T_space[1],
                    'V':rotated_K,
                   'V-I':T_space[0],
                    'Ierr':galaxy['J_err'],
                   'Verr':galaxy['K_err'],
                    'V-Ierr':T_space[1]-galaxy['K_err']
                   })
                
            # color cut to select only RGB stars
            cat = cat[((galaxy['J']-galaxy['K'])<1.3)&(((galaxy['J']-galaxy['K'])>.7))]
   
            # unrelated, but print average photometric error for TRGB error budget
            test = galaxy[(galaxy['J']<err_calcJ+.2)&(galaxy['J']>err_calcJ-.2)]
            print('average j photometric uncert:'+str(np.mean(test['J_err'])))
            print('average h photometric uncert:'+str(np.mean(test['H_err'])))
            
            cat.to_csv('/Users/abigaillee/Documents/Photometry/python_scripts/fourstar_project/fourstar.csv', index=False,\
                  sep='\t')
            
     
            # code to make a command to open the trgbpars and change values if needed
            with open("/Users/abigaillee/Documents/Photometry/python_scripts/fourstar_project/trgbpars") as f:
                lines = f.readlines()
            lines[0] = 'muguess          | '+str(mu_init)+' \n'
            with open("/Users/abigaillee/Documents/Photometry/python_scripts/fourstar_project/trgbpars", "w") as f:
                f.writelines(lines)

                
            # run command line script
            os.system('python /Users/abigaillee/Documents/Photometry/python_scripts/fourstar_project/fit_trgb.py')
            
            # now write trgb value to here so i can transform it back
            with open('/Users/abigaillee/Documents/Fourstar project/Analysis/TRGB info/TRGB_T/trgb.txt') as file:
                VALUE = float(file.read())
                
            
#             # write for t-band rectification plot in the future
#             with open('/Users/abigaillee/Documents/Fourstar project/Analysis/TRGB info/TRGB_T/'+str(name)+'_trgb_tj.txt', 'w') as f:
#                 f.write(str(VALUE))
            
#             # transform back
            x_2 = np.vstack((x,VALUE*ones,ones)) #get from TRGB code
            reverse_rotation = np.linalg.inv(A)
            final_line_K = reverse_rotation.dot(x_2) # get 'zero point' from reverse rotated line

            value_trgb = np.polyfit(final_line_K[0], final_line_K[1], deg=1)[1]

            print('trgb j='+str(value_trgb-.85))
            
            # save edge response and luminosity function
#             pd.read_csv('trgb_info.csv').to_csv('/Users/abigaillee/Photometry/python_scripts/fourstar_project/'+str(name)+'/t_j_trgblf.csv')

#             pd.read_csv('trgb_edge_info.csv').to_csv('/Users/abigaillee/Photometry/python_scripts/fourstar_project/'+str(name)+'/t_j_trgbedr.csv')

        elif band=='H':
            x = np.arange(-1,1.6,.01)
            y_int = y_intercep_H+1.62; x_int = 0
            ones = np.ones_like(x)

            # rotate data around line of slope.85
            theta = np.arctan(1.62)
            
            A = np.array([[np.cos(theta), np.sin(theta),0],[-np.sin(theta), np.cos(theta),0],[x_int-np.cos(theta)*x_int + y_int*np.sin(theta),\
                                y_int-x_int*np.sin(theta)-y_int*np.cos(theta),1]]).T
            
            ones_JH = np.ones(len(galaxy))
            
            
            JH_space = np.vstack((galaxy['J']-galaxy['K'], galaxy['H'], ones_JK))
            T_space = A.dot(JH_space)
            
            rotated_H = T_space[1] - T_space[0]
            
            # save to file for T-band TRGB plots if wanted
#             pd.DataFrame({'color':galaxy['J']-galaxy['K'], 'mag':T_space[1]}).to_csv('/Users/abigaillee/Documents/Fourstar project/Analysis/TRGB info/'+str(name)+'_tjtrgb.csv')
            
                   
            
            cat = pd.DataFrame({'I':T_space[1],
                    'V':rotated_H,
                   'V-I':T_space[0],
                    'Ierr':galaxy['H_err'],
                   'Verr':galaxy['K_err'],
                    'V-Ierr':T_space[1]-galaxy['K_err']
                   })
                
            # color cut to select only RGB stars
            cat = cat[((galaxy['J']-galaxy['K'])<1.3)&(((galaxy['J']-galaxy['K'])>.7))]
   
            
            cat.to_csv('/Users/abigaillee/Documents/Photometry/python_scripts/fourstar_project/fourstar.csv', index=False,\
                  sep='\t')
            
     
            # code to make a command to open the trgbpars and change values if needed
            with open("/Users/abigaillee/Documents/Photometry/python_scripts/fourstar_project/trgbpars") as f:
                lines = f.readlines()
            lines[0] = 'muguess          | '+str(mu_init)+' \n'
            with open("/Users/abigaillee/Documents/Photometry/python_scripts/fourstar_project/trgbpars", "w") as f:
                f.writelines(lines)

                
            # run command line script
            os.system('python /Users/abigaillee/Documents/Photometry/python_scripts/fourstar_project/fit_trgb.py')
            
            # now write trgb value to here so i can transform it back
            with open('/Users/abigaillee/Documents/Fourstar project/Analysis/TRGB info/TRGB_T/trgb.txt') as file:
                VALUE = float(file.read())
                
            
#             # write for t-band rectification plot in the future
#             with open('/Users/abigaillee/Documents/Fourstar project/Analysis/TRGB info/TRGB_T/'+str(name)+'_trgb_tj.txt', 'w') as f:
#                 f.write(str(VALUE))
            
#             # transform back
            x_2 = np.vstack((x,VALUE*ones,ones)) #get from TRGB code
            reverse_rotation = np.linalg.inv(A)
            final_line_K = reverse_rotation.dot(x_2) # get 'zero point' from reverse rotated line

            value_trgb = np.polyfit(final_line_K[0], final_line_K[1], deg=1)[1]

            print('trgb h='+str(value_trgb-1.62))
            
            # save edge response and luminosity function
#             pd.read_csv('trgb_info.csv').to_csv('/Users/abigaillee/Photometry/python_scripts/fourstar_project/'+str(name)+'/t_j_trgblf.csv')

#             pd.read_csv('trgb_edge_info.csv').to_csv('/Users/abigaillee/Photometry/python_scripts/fourstar_project/'+str(name)+'/t_j_trgbedr.csv')
            
            
        
