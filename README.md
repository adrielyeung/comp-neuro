# comp-neuro
(Course taken in spring 2019)

This course aims to understand processes for transmitting signals between neurons in the brain which leads to decision making, memory and learning, via different numerical models simulated by computers.

Code written in MATLAB. Exercises:

**Single neuron model**
- PS4: Coding the integrate-and-fire model using the Euler method to approximate, with different currents *I*.
**Neural networks**
- PS5: Coding the activity of a ring network which (with connections) provides positive feedback to close neurons and negative ones to neurons further away. Then providing an input which also depends on orientation of each neuron. Finally test how a change in input orientation affects neuron activity with and without connections.
- PS6 (majority of code provided): Simulating activity of neurons in a neural network with different mean and standard deviation in the firing rate. The relationship between synaptic weight and number of neurons is derived from the ideal behaviour in mean (low) and SD (high).
**Learning and Memory**
- PS7 (a): Implementation of the Bienenstock-Cooper-Munro (BCM) Learning Rule given two input patterns chosen randomly at the start of each simulation.
- PS7 (b): Implementation of Spike-Timing Dependent Plasticity (STDP) between two neurons, where the weight change depends on time difference from the presynaptic neuron spike.
- PS8: Implementation of the perceptron to learn patterns based on weight changes.
- PS9 (majority of code provided): Implementation of temporal difference (TD) learning, which simulates learning based on difference between predicted reward (due to a stimulus) and the actual reward.

## Acknowledgements
Credits to Dr. Claudia Clopath for teaching the course and providing example code.
