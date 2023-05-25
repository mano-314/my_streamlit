import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation 
from matplotlib.animation import FuncAnimation
from sklearn.metrics.pairwise import euclidean_distances
import sys


class initiallize() :
    def __init__(self, nat, dim ) :
        self.dim = dim
        self.nat = nat

    def vel_init(self):
        x_ = np.random.rand(self.nat, self.dim)
        y_ = np.random.rand(self.nat, self.dim)
        vel  = np.sqrt(-2*np.log(x_)) * np.cos(2 * np.pi * y_)
        return vel

    def vel_nd_lj(self):
        v = np.zeros((self.nat, self.dim))
        for i in range(self.nat) :
            a=np.random.rand(self.nat)*np.pi*2.0
            b=(np.random.rand(self.nat)-0.5)*np.pi
            if self.dim == 2 :
                v[:,0]=np.cos(b)*np.sin(a)
                v[:,1]=np.cos(b)*np.cos(a)
            elif self.dim == 3 :
                v[:,0]=np.cos(b)*np.sin(a)
                v[:,1]=np.cos(b)*np.cos(a)
                v[:,2]=np.sin(b)
        return v

    def vel_static(self):
        v = np.zeros((self.nat, self.dim))
        return v

    def pos_init(self):
        pos = np.random.rand(self.nat, self.dim)
        return pos


class potential():
    def __init__(self, pos, vel, dim, nat, dt, m ):
        self.pos = pos
        self.vel = vel
        self.dim = dim
        self.nat = nat 
        self.dt  = dt 
        self.m   = m 
        
    def print_pos(self):
        print(self.pos)

    def print_vel(self):
        print(self.vel)

    def potential_zero(self):
        pos = self.pos + self.dt*self.vel
        eu = 0
        return pos, self.vel, eu

    def lennard_jones(self, eps, sigma):
        forces = np.zeros((self.nat, self.nat))
        accel  = np.zeros((self.nat, self.dim))
        potential_u = np.zeros((self.nat, self.nat))
        for i in range(0, self.nat):
            r_ij_tmp = self.pos - self.pos[i]
            r_dis  = np.linalg.norm(r_ij_tmp , axis = 1)
            r2    = r_dis ** 2
            r_ij  = r_ij_tmp / np.linalg.norm(r_ij_tmp)

            for j_i in range(i + 1, i + self.nat):
                if j_i >= self.nat:
                    j = j_i - self.nat
                else :
                    j = j_i

                s_r_6 = ( sigma / r_dis[j] ) ** 6
                forces[i][j] = 48 * eps / r2[j] * ( s_r_6 - 0.5 ) * s_r_6
                potential_u[i][j] = 4 * eps * s_r_6 * ( s_r_6 - 1 )
                for k in range(self.dim):
                    accel[i][k] += forces[i][j] * r_ij[j][k] / self.m[i]
        self.vel += accel * self.dt
        self.pos += self.vel * self.dt + 0.5 * accel * self.dt * self.dt
        eu = np.sum(potential_u)
        return self.pos , self.vel, eu

    def coulomb(self, eps_0,nat):
        forces = np.zeros((nat, nat))
        accel  = np.zeros((nat, self.dim))
        potential_u = np.zeros((nat, nat))
        for i in range(0, nat):
            r_ij_tmp = self.pos - self.pos[i]
            r_dis  = np.linalg.norm(r_ij_tmp , axis = 1)
            r2    = r_dis ** 2
            r_ij  = r_ij_tmp / np.linalg.norm(r_ij_tmp)

            for j_i in range(i + 1, i + nat):
                if j_i >= nat:
                    j = j_i - nat
                else :
                    j = j_i

                forces[i][j] = - 0.25 /( np.pi *eps_0 * r2[j] )
                potential_u[i][j] = 0.25 /( np.pi *eps_0 * r_dis[j] )
                for k in range(self.dim):
                    accel[i][k] += forces[i][j] * r_ij[j][k] / self.m[i]
        self.vel += accel * self.dt
        self.pos += self.vel * self.dt + 0.5 * accel * self.dt * self.dt
        eu = np.sum(potential_u)
        return self.pos , self.vel, eu


    def morse(self, a_e, D_e, r_e):
        forces = np.zeros((self.nat, self.nat))
        accel  = np.zeros((self.nat, self.dim))
        potential_u = np.zeros((self.nat, self.nat))
        for i in range(0, self.nat):
            r_ij_tmp = self.pos - self.pos[i]
            r_dis  = np.linalg.norm(r_ij_tmp , axis = 1)
            r2    = r_dis ** 2
            r_ij  = r_ij_tmp / np.linalg.norm(r_ij_tmp)

            for j_i in range(i + 1, i + self.nat):
                if j_i >= self.nat:
                    j = j_i - self.nat
                else :
                    j = j_i

                forces[i][j] = 2*a_e*D_e* ( 1 -
                              np.exp( -a_e * (r_dis[j] - r_e) )) * (
                              np.exp( -a_e * (r_dis[j] - r_e) ))
                potential_u[i][j] = D_e * ( 1 -
                              np.exp( -a_e * (r_dis[j] - r_e) )) ** 2
                for k in range(self.dim):
                    accel[i][k] += forces[i][j] * r_ij[j][k] / self.m[i]
        self.vel += accel * self.dt
        self.pos += self.vel * self.dt + 0.5 * accel * self.dt * self.dt
        eu = np.sum(potential_u)
        return self.pos , self.vel, eu

    def OPPm(self, a_o, b_o, m_o, k_o, phi_o):
        forces = np.zeros((self.nat, self.nat))
        accel  = np.zeros((self.nat, self.dim))
        potential_u = np.zeros((self.nat, self.nat))
        for i in range(0, self.nat):
            r_ij_tmp = self.pos - self.pos[i]
            r_dis  = np.linalg.norm(r_ij_tmp , axis = 1)
            r2    = r_dis ** 2
            r_ij  = r_ij_tmp / np.linalg.norm(r_ij_tmp)

            for j_i in range(i + 1, i + self.nat):
                if j_i >= self.nat:
                    j = j_i - self.nat
                else :
                    j = j_i
                F_1 = -15/np.power(r_dis[j],16)
                F_2 = a_o*np.exp( -(r_dis[j] / b_o)**m_o )
                F_3 = -(m_o/b_o)*(r_dis[j]/b_o)**(m_o-1) * np.cos( k_o*r_dis[j] - phi_o )
                F_4 = k_o*np.sin( k_o*r_dis[j] - phi_o )
                forces[i][j] = F_1 + F_2 * ( F_3 + F_4 )
                potential_u[i][j] = 0
                for k in range(self.dim):
                    accel[i][k] += forces[i][j] * r_ij[j][k] / self.m[i]
        self.vel += accel * self.dt
        self.pos += self.vel * self.dt + 0.5 * accel * self.dt * self.dt
        eu = np.sum(potential_u)
        return self.pos , self.vel, eu


class periodic_condition() :

    def __init__(self,pos_in, pos_out, vel, nat, dim, dt):
        self.pos_in  = pos_in
        self.pos_out = pos_out
        self.vel     = vel
        self.nat     = nat
        self.dim     = dim

    def print_nat(self):
        print(self.nat)

    def print_dim(self):
        print(self.dim)

    def pbc_on(self):
        op5 = np.ones((self.nat,self.dim))*0.5
        pos_ = self.pos_out - np.round(self.pos_out - op5)
        return pos_ , self.vel

    def pbc_off(self):
        op5 = np.ones((self.nat,self.dim))*0.5
        y = np.abs(self.pos_out - op5) - op5
        out_in_bl  = y > 0

        out_in_int = out_in_bl#.astype(np.int64) * (-2.0) + 1.0
        vel_ = self.vel * out_in_int

        pos_ = self.pos_in + self.vel * out_in_int * self.dt
        return pos_ , vel_


def simulation(nat, infect_dist, max_step, heal_speed, dt, spread_treshold, pot_tag, outfile_name) : 

    #------- main parameters

    dim      = 2
    pbc_tag  = 1        # 0 : without pbc
                        # 1 : with pbc
    init_vel = 1        # 0 : start from static
                        # 1 : start from random velocities (Maxwell dist.)

    #------- parameters for LJ potential
    sigma = 1.0e-11
    eps   = 8
    #wm = 1.0e-9

    #------- parameters for Coulomb potential
    eps_0  = 1.0e16
    #wm     = 1.0e-5
    #------- parameters for Morse potential
    a_e    = 1.0e-3
    D_e    = 1.0e2
    r_e    = 5.0e-2
    wm     = 5.0e-5
    #------- parameters for OPPm potential
    a_o    = 0.5
    b_o    = 1.45
    m_o    = 2
    k_o    = 14.4
    phi_o  = 17.125
    # wm     = 1.0e19

    #------- control mass under potential
    #m = np.ones((nat))*wm
    m = np.random.rand(nat)*wm
    m[1] = m[1] #*10000
    # col = np.zeros((nat))
    # col[1] = 1

    #--------------------------------------
    aa = np.pi * infect_dist * infect_dist

    print( 'one circle    : ' , str( aa ) )
    print( 'area per dot  : ' , float( 1/nat ) )
    print( 'freedom       : ' , str( 1 / ( nat * aa )) )
    print( '[ def. of freedom >> area_per_dot / area_of_one_circle] ' )

    #------- CLASS


    #------ initial setting

    INIT = initiallize(nat, dim)
    pos_ = INIT.pos_init()
    if init_vel == 0 :
        vel_ = INIT.vel_static()
    elif init_vel == 1 :
        vel_ = INIT.vel_init()
    else :
        print('please set init_vel = 1 or 2')
        sys.exit()
    #vel_ = test.vel_static()

    pos_in = pos_
    vel_in = vel_
    step = 0

    ims = []
    save_infect_num = []
    healed_num = []
    status = np.ones((nat))
    experience = np.zeros((nat))


    #---------------------- graph preparation

    fig = plt.figure(figsize = (6,8))
    ax1 = fig.add_subplot(2, 1, 2)
    ax2 = fig.add_subplot(2, 1, 1)
    ax1.set_position([0.1,0.1,0.8,0.6])
    ax2.set_position([0.1,0.75,0.8,0.10])
    ax1.set_xticks([])
    ax1.set_yticks([])
    ax2.set_ylabel('number')
    ax2.set_xticks([])
    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 1)
    ax2.set_xlabel('step')

    explaination = "Number of dots : " + str(nat) + "\n" + "critical radius : " + str(infect_dist) +"\n"  +  "dt : " + str(dt) + "\n"
    explaination2 = "Recover speed : " + str(heal_speed) + "\n" + "Spread threshold : " + str(spread_treshold)

    fig.text(0.1, 0.88 , explaination + explaination2)


    #---------------------- main LOOP

    while step < max_step :

        if step == 1 :
            status[0] = 0

    #------ update positions and velocities
        pot = potential(pos_in, vel_in, dim, nat, dt, m )
        if pot_tag == 0 :
            pos_out ,vel_out, eu = pot.potential_zero(dt)
        elif pot_tag == 1 :
            pos_out ,vel_out, eu = pot.lennard_jones(eps, sigma)
        elif pot_tag == 2 :
            pos_out ,vel_out, eu = pot.coulomb(eps_0,nat)
        elif pot_tag == 3 :
            pos_out ,vel_out, eu = pot.morse(a_e, D_e, r_e)
        elif pot_tag == 4 :
            pos_out ,vel_out, eu = pot.OPPm(a_o, b_o, m_o, k_o, phi_o)
        else :
            print('pot = 0 for non-potential and 1 for LJ-potential')
            break

    #------ undate pbc
        PBC = periodic_condition(pos_in, pos_out, vel_out, nat, dim, dt)
        if pbc_tag == 0 :
            pos_out, vel_out = PBC.pbc_off()
        elif pbc_tag == 1 :
            pos_out, vel_out = PBC.pbc_on()
        else :
            print('pot = 0 for without pbc and 1 for with pbc')
            break

    #--- update status

    #------- get index for (status<spread_treshold)
        index_spread   = np.argwhere(status<spread_treshold)
    #------- get index for no experience
        index_no_exp   = np.argwhere(status==1.0)
    #------- get index for once experience >> set component in experience to 1
        index_once_exp = np.argwhere(status==0.0)
    #------- heal matrix
        heal_matrix_index = status < 1
        heal_matrix       = heal_matrix_index #.astype(np.int64)

    #-- all component which satisfy (status<spread_treshold) will be used
        for i in index_spread :
            dis_matrix = euclidean_distances(pos_out, pos_out[i])
            dis_tmp = dis_matrix > infect_dist
            dis_tmp_ = dis_tmp #.astype(np.int64)
    #------ not infected from thierself
            dis_tmp0 = dis_matrix == 0
            dis_tmp0_ = dis_tmp0 #.astype(np.int64)
    #------ only no-experience samples were considered
            status[index_no_exp] = ( status[index_no_exp]*
                                ( dis_tmp_.reshape((nat))[index_no_exp]
                                + dis_tmp0_.reshape((nat))[index_no_exp]) )

    #--- save experience : whose status once reach zeros
        experience[index_once_exp] = 1

    #--- healing with heal_speed
        status_ = status + heal_matrix * heal_speed

    #--- save number of [un-healed], and [healed]
        sum_tmp = np.sum(status_<1)
        save_infect_num.append(sum_tmp)
        healed_num.append( np.sum(experience * status_ >= 1.0) )


    #--- Visualization
        size_ = 160 / np.sqrt(nat)
        size__ = 10000 * infect_dist

        plot_circle_x = []
        plot_circle_y = []

        for index in index_spread :
            plot_circle_x.append( pos_out[index[0]][0] )
            plot_circle_y.append( pos_out[index[0]][1] )
        col = status_
    #    col[1] = 0
    #--- plot moving carier (status factor less than spread_treshold)
        im3 = ax1.scatter(plot_circle_x , plot_circle_y ,
                        s=size__, alpha=0.1 )
                        #   alpha=0.1 ,
                        #   linewidths="1.5",)
                        #   edgecolor="black" )

    #--- plot all moving dots
        im1 = ax1.scatter(pos_out[:,0], pos_out[:,1],
                        s=size_,
                        vmin=0, vmax=1.0,
                        cmap=plt.cm.viridis_r, c=np.real(col) ,
                        alpha=0.95)

    #--- plot infect number (status factor == 1)
        xx = np.linspace(0, step+1, step+1)
        ax2.set_ylim(0,nat+nat*0.1)
        col_0 = np.zeros((step+1))
        col_1 = np.ones((step+1))
        im2 = ax2.scatter(xx, save_infect_num, s=1,
                        vmin=0, vmax=1.0,
                        cmap=plt.cm.viridis_r, c=np.real(col_0) )

    #--- plot healed number (experience * status_ >= 1.0)
        im4 = ax2.scatter(xx, healed_num , s = 0.5,
                        vmin=0, vmax=1.0,
                        cmap=plt.cm.viridis_r, c=np.real(col_1) )

        ims.append([im1, im2, im3,im4])
    #   
        pos_in = pos_out
        vel_in = vel_out
        step += 1
        status = status_

    # fig.canvas.draw()
    ani = animation.ArtistAnimation(fig, ims, interval = 80 )


    fig = outfile_name + '.gif'
    ani.save( fig , writer="imagemagick")
    plt.show()
