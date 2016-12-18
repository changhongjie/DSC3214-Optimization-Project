# DSC3214 Project
# Modeling bus demand/supply
using JuMP, Gurobi

#Read in data from Project.csv
csvboard = readcsv("A1_board.csv", header = true)
boarddata = csvboard[1]
boardheader = csvboard[2]

csvalight = readcsv("A1_alight.csv", header = true)
alightdata = csvalight[1] # horizontal vector
alightheader = csvalight[2]

csvdata = readcsv("A1_data_new.csv", header = true)
data = csvdata[1]
dataheader = csvdata[2]

#Initialise data
times = boarddata[1:size(boarddata,1),1]    # exclude blank square

demand = map(Int64, boarddata[:,2:size(boarddata,2)])
alight = map(Float64, alightdata[:,2:size(alightdata,2)])

rows = size(demand,1)    # num of rows in data
columns = size(demand,2) # num of columns in data

#Parameters
capacity = map(Int64, data[1,1])                # capacity per bus
cost = map(Int64, data[1,2])                    # operating cost per bus per trip
satisfaction = map(Float64, data[1,3])          # minimum satisfaction level of customers
journeytime = map(Int64, data[1,4])             # total duration of 1 bus trip
intervaltime = map(Int64, data[1,5])            # duration of interval between despatching buses
costofbus = map(Float64, data[1,6])             # capital cost of acquiring a new bus

alighting = alightdata
journey = ceil(Int, journeytime/intervaltime)     # how many intervals is equivalent to 1 journey
M = 10000                                         # Arbitrary large number for if-else constraints
# Generate model
m = Model(solver = GurobiSolver())

@variable(m, x[1:rows] >= 0, Int)                           # decision variable on number of buses to send out
@variable(m, unsatisfied[i=1:rows, j=1:columns] <= 0, Int)  # number of customers that were unable to board bus at each stop/time
@variable(m, space[i=1:rows, j=1:(columns+1)] >= 0)         # total amount of space left on the bus(es) for the stop Sij
@variable(m, tdemand[i=1:rows, j=1:columns] >= 0, Int)      # total demand matrix which includes people who couldn't board the previous buses
@variable(m, q[1:rows, 1:columns], Bin)                     # Binary variable for if-else constraints
@variable(m, fleet >= 0, Int)                               # variable fleet size for new companies

#helper functions
function ubSpace(i,j)
  return (x[i] * capacity - space[i,j])*alighting[j] + space[i,j] - tdemand[i,j] + M*(1-q[i,j]) #set upper bound on space variable
end
function lbSpace(i,j)
  return (x[i] * capacity - space[i,j])*alighting[j] + space[i,j] - tdemand[i,j] - M*(1-q[i,j]) #set lower bound on space variable
end
function ubUnsatisfied(i,j)
  return (x[i] * capacity - space[i,j])*alighting[j] + space[i,j] - tdemand[i,j] + M*q[i,j] #set upper bound on unsatisfied variable
end
function lbUnsatisfied(i,j)
  return (x[i] * capacity - space[i,j])*alighting[j] + space[i,j] - tdemand[i,j] - M*q[i,j] #set lower bound on unsatisfied variable
end

#supply - demand constraint
for i in 1:rows
  for j in 1:columns
    #constraints on space left
    @constraint(m, space[i,j+1] <= ubSpace(i,j))
    @constraint(m, space[i,j+1] <= M*q[i,j])
    @constraint(m, space[i,j+1] >= lbSpace(i,j))
    @constraint(m, space[i,j+1] >= -M*q[i,j])
    #constraints on unsatisfied
    @constraint(m, unsatisfied[i,j] <= ubUnsatisfied(i,j))
    @constraint(m, unsatisfied[i,j] <= M*(1-q[i,j]))
    @constraint(m, unsatisfied[i,j] >= lbUnsatisfied(i,j))
    @constraint(m, unsatisfied[i,j] >= -M*(1-q[i,j]))
    #Mq constraint
    @constraint(m, (x[i] * capacity - space[i,j])*alighting[j] + space[i,j] - tdemand[i,j] <= M*q[i,j])
  end
end
# bus capacity constraint
for i in 1:rows
  @constraint(m, space[i,1] == x[i] * capacity) #nobody on the bus before 1st stop
  for j in 2:columns+1
    @constraint(m, space[i,j] <= x[i] * capacity)
  end
end
# satisfaction level constraint
@constraint(m, (-sum(unsatisfied)) <= (1-satisfaction)*sum(tdemand))
# fleet size constraint
for i in 1:(rows-journey+1)
  @constraint(m, sum{x[i+j], j=0:(journey-1)} <= fleet) #each set of j intervals meets fleet constraint
end
# move unsatisfied customers to next intervals
for j in 1:columns
  @constraint(m, tdemand[1,j] == demand[1,j]) #nobody unsatisfied before 1st interval
  for i in 1:rows-1
    @constraint(m, tdemand[i+1,j] == demand[i+1,j] - unsatisfied[i,j])
  end
end

@objective(m, Min, sum(x)*cost + fleet*costofbus) #Operating cost and Capital cost of acquiring fleet

#print(m)
solve(m)

println("Optimal Cost: ", getobjectivevalue(m))
println("Fleet size: ", getvalue(fleet))
println("Satisfaction level: ", 1-getvalue(sum(-unsatisfied))/sum(getvalue(tdemand)))
for i in 1:rows
  println("$(round(times[i])) : $(getvalue(x[i])) buses")
end

#= for debugging
for i in 1:rows
  for j in 1:columns
    println("u: $(getvalue(unsatisfied[i,j])), s: $(getvalue(space[i,j]))")
  end
end
=#
