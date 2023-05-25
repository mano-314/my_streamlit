# The-spread-of-virus


parameters to set 
- `max_step` : step for interation  
- `dt` : positions are simply updated by new_position = old_position + dt*velocity 
- `nat` : number of dots (people in community)
- `dim` : dimension of movement (basically for 2 dimension)
- `infect_dist` : critical radius for each person. 
- `heal_speed` : recover speed for sick people. 
- `spread_treshold` : how many people can transmit virus. 

All of people (dots) in this simulation were label by `status` run from 0 to 1, 
when 0 indicates completely infected status and 1 is completely healthy status. 

The incubation period is set to 0 here. It means that once healthy person gets into 
the circle created by a person who was labeled with status 
0 - `spread_treshold` (completely sick - almost sick), 
that person will be immediately infected and have one's own symptoms developed. 

The status of that person will be set to 0 immediately 
and gradually increase at the rate of `heal_speed` during the interation 
until it reach 1 (healthy status) again. Here, we assume that 
recovered people will not be infected again. 

example_200.gif was simulated by 200 dots(people). 
![](example_200.gif)

example.gif was simulated by 2000 dots(people).
![](example.gif)
