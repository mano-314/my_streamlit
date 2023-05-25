import streamlit as st 
import numpy as np 
from simulation_corona import simulation
import os 

st.markdown(f'<b style="color:black;font-size:25px;">Does social distancing work?</b>', unsafe_allow_html=True)

col1, col2, col3, col4 = st.columns(4)

with col1:
    nat = st.slider(" " , 0, 200 , 100 )
    st.metric("Number of people", nat) 
    st.markdown("Number of people in society in modeling")

with col2:
    max_step=st.slider(" " ,0,500 , 100 )
    st.metric("Number of steps", max_step ) 
    st.markdown("How long the simulation?")

with col3:
    dt=st.slider(" " ,1,100, 50 )
    st.metric("Time step", dt ) 
    st.markdown("Time step for simulation > The speed of people's moving")
    dt = dt*0.0001

with col4:
    pot=st.selectbox(" ", ["non-potential","LJ","Coulomb"],index=1)
    st.metric("Select interaction", pot ) 
    st.markdown("People always have relationship, sometimes get closer, sometimes takes distance.")


col1, col2, col3, col4 = st.columns(4)

with col1:
    infect_dist=st.slider(" ", 0.00,0.10, 0.06 , step = 0.01 ) 
    st.metric("Infect distance", infect_dist ) 
    st.markdown("If you get inside the radius of someone who got infected, you will be the next one")

with col2: 
    heal_speed=st.slider(" ",0.0,0.1, 0.01)
    st.metric("heal speed", heal_speed ) 
    st.markdown("The recovering speed after get infected (assume that those who experienced once will not get infected again)")

with col3:
    spread_treshold=st.slider(" " ,0.0,1.0, 0.5)
    st.metric("Spread treshold", spread_treshold ) 
    st.markdown("How many people can transmit virus during got infected?")
    
    
with col4:
    outfile_name=st.text_input(" " , "output")
    st.metric("Output file name", outfile_name )  
    st.markdown("Output filename for the simulation")

    if pot=="non-potential":
        pot_tag=0
    elif pot=="LJ":
        pot_tag=1
    elif pot=="Coulomb":
        pot_tag=2
    elif pot=="Morse":
        pot_tag=3
    else:
        pass 

st.write('\n')
st.write('\n')

col1, col2, col3, col4, col5 = st.columns(5)
with col2:
    st.write('\n')
    lsim = st.button("RUN simulation")
with col4:
    example = st.selectbox("Select example" , ["small", "medium","large","output"] , index = 0)
    
st.write('\n')
st.write('\n')



if lsim:
    process = st.text("processing ... \nprocessing ... \nprocessing ...")
    simulation(nat, infect_dist, max_step, heal_speed, dt, spread_treshold, pot_tag, outfile_name)
    
    process.text("DONE")
    fig = outfile_name + '.gif'
    
    st.image(fig,width = 700)

    with open(fig, "rb") as file:
        btn = st.download_button(
            label="Download",
            data=file,
            file_name=fig,
            mime="image/gif"
        )

else : 
    if example=="output": 
        fig = outfile_name + '.gif'
    else :
        fig = 'example/example_{}.gif'.format(example)
    if os.path.exists(fig):

        st.image(fig,width = 700)
        with open(fig, "rb") as file:
            btn = st.download_button(
                label="Download",
                data=file, 
                file_name=fig,
                mime="image/gif"
            )
    