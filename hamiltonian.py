#1:	Start finding solution at 1
#2: Try to quickly find the next solution based on the last one
#3: If that fails, try finding the path using a brute force method
#4:	After a solution is found (on step 2 or 3) print it, and increment n, goto #1
#5:	Otherwise, print why it failed and then goto #1.

def main():
	curentSize = 0
	solList = []
	lastSolution = None
	solFound = False
	while True:
		print("Finding solution for",curentSize+1)
		if(lastSolution is not None):
			solFound = False
			#Generate Square Graph to be used in find next
			squareGraph = generateSquareGraph(curentSize+1)
			for s in getNextPerm(lastSolution):
				n = findNext(s,squareGraph)
				if(n is not None):
					if(not sanityCheck(n)):
						print("Error, failing sanity check")
					print("Next solution found for",curentSize+1)
					solList.append([curentSize+1,n])
					lastSolution = n
					curentSize += 1
					solFound = True
					print(lastSolution)
					input()
					break
			if(solFound):
				continue
			else:
				print("No quick solution for",curentSize+1,"was found")
		print("starting brute force method")
		lastSolution = hamiltonianPath(generateSquareGraph(curentSize+1),lastSolution)[0]
		solList.append([curentSize+1,lastSolution])
		print(lastSolution)
		if(type(lastSolution) == type("")):
			lastSolution = None
		curentSize += 1
		input()

#Creates a graph for all numbers up to n, which relates two number if their sum
#	is a square number. The edges are bidirectional, and represented on both numbers.
def generateSquareGraph(n):
	#Generate List of Squares less than n
	squares = []
	i = 1
	while(i**2 <= 2*n):
		squares.append(i**2)
		i += 1
	
	graphDict = {}
	
	for i in range(1,n+1):
		l = []
		for s in squares:
			if(s > i and s-i <= n and s-i != i):
				l.append(s-i)
		graphDict[i] = set(l)
	return graphDict	

#Ensures a solution is valid
def sanityCheck(sol):
	n = max(sol)
	
	#Generate the relevant list of squares
	squares = []
	i = 1
	while(i**2 <= 2*n):
		squares.append(i**2)
		i += 1
	squares = set(squares)
	
	#Check all items add to a square
	for i in range(len(sol)-1):
		if((sol[i] + sol[i+1]) not in squares):
			return False
			
	#Check all numbers are represented exactly once
	counts = [0]*n
	for i in sol:
		if(i < 0 or i > n):
			return False
		counts[i-1] += 1
	for c in counts:
		if c != 1:
			return False
	
	return True

#Takes string representation of a solution and returns the list
def toList(sol):
	sol = sol.replace("[","")
	sol = sol.replace("]","")
	return [int(x) for x in sol.split(", ")]

#Takes list representation of a solution and returns the hashable string
def toString(sol):
	s2 = [str(x) for x in sol]
	return ", ".join(s2)

#Generates more solutions from a single solution by performing
#	mutating the list by pairing and end of the list to a matching
#	element on the inside. The list is cut at the match, one portion is reversed
#	and the new lists are stitched together by the matching element and the end it matched.

#This is useful for creating new solutions for the "findNext" solution to attempt
#	to quickly solve. It turns out, as the size of the list increases, the odds
#	this method of finding the next solution works should DRAMATICALLY increase.
def getNextPerm(sol):
	frontier = set([toString(sol)])
	visited = set()
	n = max(sol)
	g = generateSquareGraph(n)
	
	squares = []
	i = 1
	while(i**2 <= 2*n):
		squares.append(i**2)
		i += 1
	
	newSolutions = []
	while len(frontier) > 0:
		nextStr = frontier.pop()
		visited.add(nextStr)
		next = toList(nextStr)
		
		first = next[0]
		last = next[-1]
		
		firstMatch = g[first]
		lastMatch = g[last]
		
		for f in firstMatch:
			after = next.index(f)
			a = next[after:]
			a.reverse()
			newSolutions.append(a + next[:after])
		
		for f in lastMatch:
			before = next.index(f)+1
			a = next[:before]
			a.reverse()
			newSolutions.append(next[before:] + a)
		
		for s in newSolutions:
			stringS = toString(s)
			if(stringS not in visited and stringS not in frontier):
				frontier.add(stringS)
				yield s
				
	raise StopIteration

#Sometimes you can just stick n+1 to the end of a solution for n and your done.
#	This is more effective when you can also generate 
#	new solutions a single one, and then check it with this.
#
#Because of how this method works, if a certain list for n has no solutions where n is
#	on either the first or last element, then this will fail to find a solution.
#	I'm fairly sure this only happens on 27,28, and 29. After that happens, this method
#	appears to always work
#
#It also looks 1 cut, reverse, and stitch ahead (method used in genPerm).
#	This probably speeds things up maybe.
def findNext(lastSolList,squareGraph):
	n = max(lastSolList) + 1
	
	nEdges = squareGraph[n]
	
	first = lastSolList[0]
	last = lastSolList[-1]
	
	squares = []
	i = 1
	while(i**2 <= 2*n):
		squares.append(i**2)
		i += 1
	
	
	#Place at start or end if possible
	if(lastSolList[0] in nEdges):
		return [n] + lastSolList
	if(lastSolList[-1] in nEdges):
		return lastSolList + [n]
		
	#Find neighbors to edges, try stitching to end
	for e in nEdges:
		eIndex = lastSolList.index(e)
		before = lastSolList[eIndex - 1]
		after = lastSolList[eIndex + 1]
		
		if((before+last) in squares):
			a = lastSolList[:eIndex]
			a.reverse()
			return [n] + lastSolList[eIndex:] + a
			
		if((after+first) in squares):
			a = lastSolList[eIndex+1:]
			a.reverse()
			return a + lastSolList[:eIndex+1] + [n]
	return None

#
#ALL FUNCTIONS BELOW THIS POINT FOR ARE BRUTE FORCE METHOD
#This is really only used between 1 and 30 when iterating solutions
#	by generating solutiosn and "findNext"ing them doesn't really work yet
	
#Function to detect if graph is well connected
#Poorly connected graphs cannot have a hamiltonian path by definition
def isStronglyConnected(graphDict):
	visited = {}
	frontier = {}
	current = None
	ret = True
	
	#Choose single node
	for k in graphDict:
		current = k
	
	visited[current] = True
	
	#Add nodes to frontier
	for v in graphDict[k]:
		frontier[v] = True
		
	while len(frontier) > 0:
		
		#Get any element from frontier
		for k,v in frontier.items():
			current = k
			break
		
		
			
		visited[current] = True
		for v in graphDict[current]:
			#Add egde to frontier if not yet added or visited
			if v not in visited and v not in frontier:
				frontier[v] = True
				
		#Remove that element from frontier
		del(frontier[current])
		if(current in frontier):
			print("ERRRROR")
	
	#Check that all nodes in the graph are also visited
	for k in graphDict:
		if k not in visited:
			ret = False
			break
			
	return ret

#Checks that no subgraph is disconnected from the main graph
#Used as an early exit condition during a hamiltonian brute force
def isStronglyConnectedHamiltonian(graphDict,visitedHamiltonian,head,tail):
	visited = {}
	frontier = {}
	if(head != tail):
		frontier[tail] = True
	current = None
	ret = True
	
	current = head
	
	visited[current] = True
	
	#Add nodes to frontier
	for v in graphDict[current]:
		if(v not in visitedHamiltonian):
			frontier[v] = True
		
	while len(frontier) > 0:
		#Get any element from frontier
		for k,v in frontier.items():
			current = k
			break
			
		visited[current] = True
		for v in graphDict[current]:
			#Add egde to frontier if not yet added or visited
			if v not in visited and v not in frontier and v not in visitedHamiltonian:
				frontier[v] = True
				
		#Remove that element from frontier
		del(frontier[current])
	
	#Check that all nodes in the graph are also visited
	for k in graphDict:
		if k not in visited and k not in visitedHamiltonian:
			ret = False
			break
			
	return ret
	

def hamiltonianBacktrack(graphDict):
	#Create ordered list from least connected nodes to most connected nodes
	vertices = []
	for k in graphDict:
		vertices.append(k)
	
	#Order the verticies from least connected to most connected
	#
	#I suspect starting with the least connected nodes will increase the speed
	#	of finding a solution, but I have no evidence for this claim
	vertices.sort(key=lambda x: len(graphDict[x]))
	
	#Identify which nodes are leaves
	leaves = []
	for k in graphDict:
		if(len(graphDict[k]) == 1):
			leaves.append(k)
	
	#If leaves exist, only try starting from leaves
	#	Valid solutions must start and end with a leaf if
	#	leaves are present
	if(len(leaves) > 0):
		for v in leaves:
			r = hamiltonianBacktrackExt(graphDict,vertices,v,v,set([v]),[v])
			if(r[1]):
				return r
	else:
		#otherwise, Try starting from every vertex instead
		for v in vertices:
			r = hamiltonianBacktrackExt(graphDict,vertices,v,v,set([v]),[v])
			if(r[1]):
				return r
	
	#Fail to solve
	return ("No path found", False)
	
	

def hamiltonianBacktrackExt(graphDict,vertices,head,tail,visited,path):
	#print(head,tail,visited)
	
	#Check if done
	if(len(visited) == len(graphDict)):
		return (path,True)
	
	#Check for graph connectivity
	if(not isStronglyConnectedHamiltonian(graphDict,visited,head,tail)):
		return ("Failed connectivity check",False)
	
	#Check if head has any pathable nodes
	for v in graphDict[head]:
		if v not in visited:
			newVisited = set(visited)
			newVisited.add(v)
			r = hamiltonianBacktrackExt(graphDict,vertices,v,tail,newVisited,path + [v])
			if(r[1]):
				return r
				
	#Check if tail has any pathable nodes
	for v in graphDict[tail]:
		if v not in visited:
			newVisited = set(visited)
			newVisited.add(v)
			r = hamiltonianBacktrackExt(graphDict,vertices,head,v,newVisited,[v] + path)
			if(r[1]):
				return r
				
	#Fail to solve
	return ("No path found",False)
	
#Wrapper function for the hamiltonian path finding functions
#Does some checks on the graph to determine solution impossibility faster
#
#The tuple return type for this is just horrible, but its used so little in the program
#	that I don't even care anymore.
def hamiltonianPath(graphDict,lastSol):
	retMessage = "No message was set"
	retSuccess = None
	
	#Check Graph Connectivity Constraint
	if(not isStronglyConnected(graphDict)):
		retMessage = "Graph is not strongly connected, no hamiltonian path exists"
		retSucess = False
		return (retMessage,retSuccess)
	
	
	#Check Leaf Constraint
	leafCount = 0
	for k,v in graphDict.items():
		if(len(v) == 1):
			leafCount += 1
	
	#Fail on more than 2 leaves
	if(leafCount > 2):
		retMessage = "Graph has " + str(leafCount) + " leaves which is more than 2, no hamiltonian path exists"
		retSuccess = False
		return (retMessage,retSuccess)
	
	#Try to solve fast with mutations if last sol is given
	if(type(lastSol) == type([])):
		m = findNext(lastSol,generateSquareGraph(max(lastSol)+1))
		if m is not None:
			return (m,True)
	
	#Perform backtracking algorithm to brute force valid paths
	r = hamiltonianBacktrack(graphDict)
	
	return r

#Python meme
if __name__ == "__main__":
	main()