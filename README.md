Main components:

RandomNormal function: This function generates random numbers from a normal distribution using the Polar Marsaglia method. It's used to simulate the stochastic nature of the oscillator.
Main function: This is the entry point of the program. It takes command-line arguments for the oscillator's parameters: alpha, sigma, M (time steps), N (iterations), and an optional seed for the random number generator.
Simulation:

The program initializes the oscillator's state with random values for X and Y.
It then iterates N times, updating the oscillator's state according to the stochastic differential equations:
X is updated based on Y and the oscillator's parameters.
Y is updated based on X, Y, and the oscillator's parameters, with a random component added using the RandomNormal function.
At each iteration, the program checks if the oscillator's state is within a certain limit (lim) of the equilibrium points (±alpha). If so, it increments a probability counter (p) at specific time points (every 100 time steps).
The program writes the probability distribution (p) to a file named Prob.out.
Output:

The program prints the input parameters and the time required to run the simulation.
It writes the probability distribution to the Prob.out file.
