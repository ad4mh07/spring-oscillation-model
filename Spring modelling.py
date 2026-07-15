import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import matplotlib.animation as animation

#Original problem


#np.sqrt((T_A+T_C)/m)

def spring_model(AC, T_A, NL_A, T_C, NL_C, m, k, offset, A=0, omega=0, g=9.8, dp=2, time=[0,10],plot=True, animate=False, energy=False):

    ex = AC - NL_A - NL_C #The sum of the extensions
    e_A = round( (T_C * ex + m * g) / (T_A + T_C) , 2 )
    e_C = ex - e_A


    a = round (k / m, 2)
    b = round( (T_A + T_C) / m , 2)
    c = round(-m*g + T_A * e_A - T_C * e_C, 2)

    if A == 0 and omega == 0:
        def f(t, A, omega):
                return A* np.sin(omega * t)

        def system(t, state, a, b):
            X, v = state
            dX_dt = v
            dv_dt = -a * v - b * X + f(t, A, omega)
            
            return [dX_dt, dv_dt]
        
    else:
         def system(t, state, a, b):
            X, v = state
            dX_dt = v
            dv_dt = -a * v - b * X
        
            return [dX_dt, dv_dt]

    t_eval = np.linspace(time[0], time[1], 100* (time[1]-time[0]) )

    initial_conditions = [offset,0]

    sol = solve_ivp(system, time, initial_conditions, t_eval=t_eval, args=(a, b))

    X_values = sol.y[0]  # This is X(t)
    v_values = sol.y[1]  # This is X'(t)


    print(f"\nAt equilibirum, the mass is at a height of {NL_C+e_C} m. \nThe spring AB is extended {e_A:.2f} m & spring BC is extended {e_C:.2f} m. \n")
    print(f"The displacement is modeeled by the equation:\ndV + {a} V + {b} X + c = {A} sin({omega}x)\n")
    if plot==True:
        print(f"The spring is pulled {offset} m downwards, resulting in this motion:")
        plt.plot(t_eval, X_values)
        plt.xlabel('Time')
        plt.ylabel('Displacement of mass below equilibium')
        plt.title("Dampened movement of a mass connected to a spring above & below.")
        plt.show()

    if animate==True:
        y=NL_C+e_C

        fig, ax=plt.subplots()
        ax.set_xlim(-1,1)
        ax.set_ylim(-0.5 , AC+0.5)
        ax.plot(0,NL_C+e_C, 'x', color='black')


        def update(frame):
            ax.clear()

            ax.set_xticks([])
            ax.set_title("Animation of the spring system")
            ax.set_ylabel('Height')

            ax.set_xlim(-1, 1)
            ax.set_ylim(-0.5 , AC+0.5)

            ax.plot([0, 0], [AC , 0], color='gray') #Vertical line, represents the springs
                    
            ax.plot([-0.3,0.3],[AC , AC], color='black') #plots a continuous line from start point to end point
            ax.plot([-0.3,0.3],[0,0], color='black') # [x1,x2],[y1,y2]

            ax.text(-0.02,AC+0.03,'A')
            ax.text(-0.02,-0.31,'C')

            ax.plot(0,NL_C+e_C, 'x', color='black')
            ax.text(-0.4,NL_C+e_C-0.1,'Equilibium')

            ax.plot(0,NL_C+e_C-offset, 'o', color='black')
            ax.text(-0.5,NL_C+e_C-offset-0.1,'Starting point')

            ax.plot(0, -X_values[frame]+NL_C+e_C, 'o', color='red')


        ani=animation.FuncAnimation(fig=fig, func=update, frames=len(X_values), interval=10)
        plt.show()  
         
    def E_k(x):
         return( 0.5*m*v_values**2 )
    
    def E_g(x):
         return( m*g*(NL_C+e_C-X_values) )
    
    def E_ep(x):
         return( 0.5*T_C*(e_C-X_values)**2 + 0.5*T_A*(e_A+X_values)**2 )

    def E_lost(x):
         dt = t_eval[1] - t_eval[0]
         return(np.cumsum(a * v_values**2) * dt)
    
    if energy==True:
        fig, ax=plt.subplots()

        plt.plot(t_eval, E_k(t_eval)+E_g(t_eval)+E_ep(t_eval)-E_lost(t_eval), color='black', label='Total energy')
        plt.plot(t_eval, E_k(t_eval), color='blue', label='Kinetic energy')
        plt.plot(t_eval, E_g(t_eval), color='Green', label='Gravitational potential energy')
        plt.plot(t_eval, E_ep(t_eval), color='orange', label='Elastic potential energy')
        plt.plot(t_eval, E_lost(t_eval), color='red', label='Dissipated energy (cumulative)')

        plt.legend(loc="upper right", fontsize = 8)
        plt.title("Energy in the system")

        ax.set_xlabel("Time")
        ax.set_ylabel("Energy (joules)")
        plt.show()
    


#spring_model(AC=8, T_A=8, NL_A=2, T_C=3, NL_C=3, m=0.5, k=0.2, A=8, omega=4.7,  offset=2, time=[0,20],plot=True, animate=True, energy=True)

spring_model(AC=8, T_A=8, NL_A=2, T_C=3, NL_C=3, m=2, k=0.5, offset=0.6, A=5, omega=4.7, animate=True, energy=True)


