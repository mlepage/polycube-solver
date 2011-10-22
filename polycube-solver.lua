-- Polycube cubic dissection solver
-- Copyright (C) 2009 Marc Lepage

require "bitmatrix"
local mcreate, mcopy, mget, mset, mprint = bitmatrix.create, bitmatrix.copy, bitmatrix.get, bitmatrix.set, bitmatrix.print
local minsertrow, mremoverow, minsertcol, mremovecol = bitmatrix.insertrow, bitmatrix.removerow, bitmatrix.insertcol, bitmatrix.removecol
local abs, ceil, rep, tremove, write = math.abs, math.ceil, string.rep, table.remove, io.write

-- Dimensions of box (width, height, depth)
local W, H, D

-- Used for handling solutions
local solline, solcount = "", 0

-- Copies a table (array values only)
local function tcopy(t)
	local tc = {}
	for i, v in ipairs(t) do
		tc[i] = v
	end
	return tc
end

-- Converts a 0-based xyz coord to a 1-based column index
local function xyz2j(x, y, z)
	-- Z-major order (X varies fastest)
	return 1 + z*H*W + y*W + x
end

-- Functions to perform the 24 axis-to-axis orientations
local orifuncs =
{
	-- +z is up, rotate around up
	function(x, y, z) return  x,  y,  z end,
	function(x, y, z) return -y,  x,  z end,
	function(x, y, z) return -x, -y,  z end,
	function(x, y, z) return  y, -x,  z end,
	-- +y is up, rotate around up
	function(x, y, z) return  x, -z,  y end,
	function(x, y, z) return  z,  x,  y end,
	function(x, y, z) return -x,  z,  y end,
	function(x, y, z) return -z, -x,  y end,
	-- +x is up, rotate around up
	function(x, y, z) return -z,  y,  x end,
	function(x, y, z) return -y, -z,  x end,
	function(x, y, z) return  z, -y,  x end,
	function(x, y, z) return  y,  z,  x end,
	-- -z is up, rotate around up
	function(x, y, z) return -x,  y, -z end,
	function(x, y, z) return -y, -x, -z end,
	function(x, y, z) return  x, -y, -z end,
	function(x, y, z) return  y,  x, -z end,
	-- -y is up, rotate around up
	function(x, y, z) return  x,  z, -y end,
	function(x, y, z) return -z,  x, -y end,
	function(x, y, z) return -x, -z, -y end,
	function(x, y, z) return  z, -x, -y end,
	-- -x is up, rotate around up
	function(x, y, z) return -z, -y, -x end,
	function(x, y, z) return  y, -z, -x end,
	function(x, y, z) return  z,  y, -x end,
	function(x, y, z) return -y,  z, -x end,
}

-- Named pieces
local pieces =
{
	-- monocube
	{ name="1_", cubes = { {0,0,0} } },
	-- dicube
	{ name="2_", cubes = { {0,0,0},{0,1,0} } },
	-- tricubes
	{ name="3I", cubes = { {0,0,0},{0,1,0},{0,2,0} } },
	{ name="3L", cubes = { {0,1,0},{0,0,0},{1,0,0} } },
	-- tetracubes (2D tetrominoes)
	{ name="4I", cubes = { {0,0,0},{0,1,0},{0,2,0},{0,3,0} } },
	{ name="4O", cubes = { {0,0,0},{0,1,0},{1,1,0},{1,0,0} } },
	{ name="4L", cubes = { {0,2,0},{0,1,0},{0,0,0},{1,0,0} } },
	{ name="4S", cubes = { {0,2,0},{0,1,0},{1,1,0},{1,0,0} } },
	{ name="4T", cubes = { {0,1,0},{1,1,0},{2,1,0},{1,0,0} } },
	-- tetracubes (3D)
	{ name="4^", cubes = { {0,1,0},{0,0,0},{1,0,0},{0,0,1} } },
	{ name="4<", cubes = { {0,1,0},{0,0,0},{1,0,0},{0,1,1} } },
	{ name="4>", cubes = { {0,1,0},{0,0,0},{1,0,0},{1,0,1} } },
	-- pentacubes (2D pentominoes)
	{ name="F_", cubes = { {2,2,0},{1,2,0},{1,1,0},{1,0,0},{0,1,0} } },
	{ name="I_", cubes = { {0,4,0},{0,3,0},{0,2,0},{0,1,0},{0,0,0} } },
	{ name="L_", cubes = { {0,3,0},{0,2,0},{0,1,0},{0,0,0},{1,0,0} } },
	{ name="N_", cubes = { {0,0,0},{0,1,0},{1,1,0},{1,2,0},{1,3,0} } },
	{ name="P_", cubes = { {0,0,0},{0,1,0},{0,2,0},{1,2,0},{1,1,0} } },
	{ name="T_", cubes = { {0,2,0},{1,2,0},{2,2,0},{1,0,0},{1,1,0} } },
	{ name="U_", cubes = { {0,1,0},{0,0,0},{1,0,0},{2,0,0},{2,1,0} } },
	{ name="V_", cubes = { {0,2,0},{0,1,0},{0,0,0},{1,0,0},{2,0,0} } },
	{ name="W_", cubes = { {0,2,0},{0,1,0},{1,1,0},{1,0,0},{2,0,0} } },
	{ name="X_", cubes = { {0,1,0},{1,1,0},{2,1,0},{1,2,0},{1,0,0} } },
	{ name="Y_", cubes = { {0,3,0},{0,2,0},{0,1,0},{0,0,0},{1,1,0} } },
	{ name="Z_", cubes = { {0,2,0},{1,2,0},{1,1,0},{1,0,0},{2,0,0} } },
	-- pentacubes (3D with mirror symmetry)
	{ name="Q_", cubes = { {0,0,0},{0,1,0},{1,1,0},{1,0,0},{1,0,1} } },
	{ name="A_", cubes = { {0,1,1},{0,1,0},{0,0,0},{1,0,0},{1,0,1} } },
	{ name="T1", cubes = { {0,1,0},{1,1,0},{2,1,0},{1,0,0},{1,1,1} } },
	{ name="T2", cubes = { {0,1,0},{1,1,0},{2,1,0},{1,0,0},{1,0,1} } },
	{ name="L3", cubes = { {0,2,0},{0,1,0},{0,0,0},{1,0,0},{0,0,1} } },
	-- pentacubes (3D chiral pairs)
	{ name="L1", cubes = { {0,2,0},{0,1,0},{0,0,0},{1,0,0},{0,2,1} } },
	{ name="J1", cubes = { {1,2,0},{1,1,0},{1,0,0},{0,0,0},{1,2,1} } },
	{ name="L2", cubes = { {0,2,0},{0,1,0},{0,0,0},{1,0,0},{0,1,1} } },
	{ name="J2", cubes = { {1,2,0},{1,1,0},{1,0,0},{0,0,0},{1,1,1} } },
	{ name="L4", cubes = { {0,2,0},{0,1,0},{0,0,0},{1,0,0},{1,0,1} } },
	{ name="J4", cubes = { {1,2,0},{1,1,0},{1,0,0},{0,0,0},{0,0,1} } },
	{ name="N1", cubes = { {1,2,0},{1,1,0},{0,1,0},{0,0,0},{1,2,1} } },
	{ name="S1", cubes = { {0,2,0},{0,1,0},{1,1,0},{1,0,0},{0,2,1} } },
	{ name="N2", cubes = { {1,2,0},{1,1,0},{0,1,0},{0,0,0},{1,1,1} } },
	{ name="S2", cubes = { {0,2,0},{0,1,0},{1,1,0},{1,0,0},{0,1,1} } },
	{ name="V1", cubes = { {1,1,0},{1,0,0},{2,0,0},{1,1,1},{0,1,1} } },
	{ name="V2", cubes = { {1,1,0},{1,0,0},{0,0,0},{1,1,1},{2,1,1} } },
}

-- Gets a piece by name
local function getpiece(name)
	for _, piece in ipairs(pieces) do
		if piece.name == name then return piece end
	end
end

-- Adds a piece to the problem matrix by calculating all positions for one or all
-- orientations and removing any non-unique placements
local function addpiece(mat, piece, constrain_x, constrain_y, constrain_z, lockcount)
	local firstrow = mat.m+1
	minsertcol(mat, mat.n+1)
	mat.hdr[#mat.hdr+1] = piece.name
	mat.count[#mat.count+1] = 0

	-- Find extents of piece
	local xmin, xmax, ymin, ymax, zmin, zmax = 999999999, 0, 999999999, 0, 999999999, 0
	for _, cube in ipairs(piece.cubes) do
		if cube[1] < xmin then xmin = cube[1] end
		if xmax < cube[1] then xmax = cube[1] end
		if cube[2] < ymin then ymin = cube[2] end
		if ymax < cube[2] then ymax = cube[2] end
		if cube[3] < zmin then zmin = cube[3] end
		if zmax < cube[3] then zmax = cube[3] end
	end

	for _, orifunc in ipairs(orifuncs) do
		-- Re-orient extents
		local xmin, ymin, zmin = orifunc(xmin, ymin, zmin)
		local xmax, ymax, zmax = orifunc(xmax, ymax, zmax)
		if xmax < xmin then xmin, xmax = xmax, xmin end
		if ymax < ymin then ymin, ymax = ymax, ymin end
		if zmax < zmin then zmin, zmax = zmax, zmin end

		-- Calculate number of positions in each direction
		local xp, yp, zp = W - (xmax-xmin), H - (ymax-ymin), D - (zmax-zmin)
		if constrain_x then xp = ceil(xp/2) end
		if constrain_y then yp = ceil(yp/2) end
		if constrain_z then zp = ceil(zp/2) end

		-- For each position
		for zo = -zmin, -zmin + zp - 1 do
			for yo = -ymin, -ymin + yp - 1 do
				for xo = -xmin, -xmin + xp - 1 do
					-- Assume unique
					minsertrow(mat, mat.m+1)
					mset(mat, mat.m, mat.n)
					for _, cube in ipairs(piece.cubes) do
						local x, y, z = orifunc(cube[1], cube[2], cube[3])
						x, y, z = xo+x, yo+y, zo+z
						assert(0 <= x and x < W)
						assert(0 <= y and y < H)
						assert(0 <= z and z < D)
						mset(mat, mat.m, xyz2j(x, y, z))
					end
					
					-- Check assumption
					local unique = true
					for i = firstrow, mat.m-1 do
						local row = mat[mat.m]
						unique = false
						for k = 1, #row do
							if mat[i][k] ~= row[k] then unique = true; break end
						end
						if not unique then mremoverow(mat, mat.m); break end
					end
					-- Update count
					if unique then
						for j = 1, mat.n do
							if mget(mat, mat.m, j) == 1 then
								mat.count[j] = mat.count[j] + 1
							end
						end
					end
				end
			end
		end
		if _ == lockcount then break end
	end
end

-- Prints a solution
local function printsol(sol)
	for y = H-1, 0, -1 do
		print(solline)
		for z = 0, D-1 do
			if z ~= 0 then write("  ") end
			for x = 0, W-1 do
				if x == 0 then write("|") end
				local name = "  "
				local j = xyz2j(x, y, z)
				for i = 1, sol.m do
					if mget(sol, i, j) == 1 then
						for j = W*H*D+1, sol.n do
							if mget(sol, i, j) == 1 then
								name = sol.hdr[j]
								break;
							end
						end
						break
					end
				end
				write(name, "|")
			end
		end
		write("\n")
		if y == 0 then print(solline) end
	end
end

-- Handles a solution
local function foundsol(sol)
	solcount = solcount + 1
	print("Solution " .. solcount)
	printsol(sol)
end

-- Knuth's Algorithm X
-- http://en.wikipedia.org/wiki/Algorithm_X
local function solve(mat, sol, solfunc)
	-- 1) If the matrix is empty, the problem is solved, terminate successfully
	if mat.n == 0 then
		solfunc(sol)
		return
	end

	-- 2) Otherwise, choose a column (deterministically)
	local minj, mincount = 0, 999999999
	for j = 1, mat.n do
		if mat.count[j] < mincount then minj, mincount = j, mat.count[j] end
	end
	local c = minj
	
	if mincount == 0 then
		-- terminate unsuccessfully
		return
	end
	
	-- 3) Choose a row such that the value at the row and column is 1 (nondeterministically)
	for i = 1, mat.m do
		if mget(mat, i, c) == 1 then
			local r = i
			
			-- Operate on copies
			local oldmat, oldsol = mat, sol
			local mat, sol = mcopy(mat), mcopy(sol)
			mat.hdr, mat.count = tcopy(oldmat.hdr), tcopy(oldmat.count)
			sol.hdr = oldsol.hdr -- OK to copy by reference
			
			-- 4) Include the chosen row in the partial solution
			minsertrow(sol, sol.m+1)
			for j = 1, mat.n do
				if mget(mat, r, j) == 1 then
					local name = mat.hdr[j]
					for j = 1, sol.n do
						if sol.hdr[j] == name then
							mset(sol, sol.m, j, 1)
							break
						end
					end
				end
			end
			
			-- 5) Reduce the matrix
			for j = oldmat.n, 1, -1 do
				if mget(oldmat, r, j) == 1 then
					for i = mat.m, 1, -1 do
						if mget(mat, i, j) == 1 then
							for j = 1, mat.n do
								if mget(mat, i, j) == 1 then
									mat.count[j] = mat.count[j] - 1
								end
							end
							mremoverow(mat, i)
						end
					end
					mremovecol(mat, j)
					tremove(mat.hdr, j)
					tremove(mat.count, j)
				end
			end
			
			-- 6) Repeat this algorithm recursively on the reduced matrix
			solve(mat, sol, solfunc)
		end
	end
end

-- Load the problem
if (not arg[1]) then print("usage: lua polycube-solver.lua <problem.txt>") os.exit(1) end
loadfile(arg[1])()
print("Polycube Solver 1.0")

-- Width, height, depth of the box to fill
W, H, D = problem.box[1], problem.box[2], problem.box[3]
print("Box dimensions: " .. W .. "x" .. H .. "x" .. D)

-- Create the problem matrix (with header and counts)
local mat = mcreate(0, W*H*D)
mat.hdr = {}
for z = 0, W-1 do
	for y = 0, H-1 do
		for x = 0, D-1 do
			mat.hdr[#mat.hdr+1] = x .. y .. z
		end
	end
end
mat.count = {}
for j = 1, mat.n do
	mat.count[#mat.count+1] = 0
end

-- Add pieces to problem matrix
write("Pieces:")
for _, name in ipairs(problem.pieces) do
	write(" " .. name)
	local piece = getpiece(name)
	local constrain_x = problem.constrain == name or problem.constrain_x == name
	local constrain_y = problem.constrain == name or problem.constrain_y == name
	local constrain_z = problem.constrain == name or problem.constrain_z == name
	addpiece(mat, piece, constrain_x, constrain_y, constrain_z, problem.lock == name and (problem.lockcount or 1) or 24)
end
print("")
if problem.contrain or problem.constrain_x or problem.constrain_y or problem.constrain_z then
	print("Constraints: " ..
		(problem.constrain or "--") .. " " ..
		(problem.constrain_x or "--") .. " " ..
		(problem.constrain_y or "--") .. " " ..
		(problem.constrain_z or "--"))
end
if problem.lock then
	print("Locked: " .. problem.lock .. " (" .. (problem.lockcount or 1) .. ")")
end

-- Print the matrix if desired
print("Problem matrix: " .. mat.m .. "x" .. mat.n)
mprint(mat)

-- Create the solution matrix (with header) and prepare it for printing
local sol = mcreate(0, mat.n)
sol.hdr = tcopy(mat.hdr)
do
	local frag = "+" .. rep("--+", W)
	solline = frag .. rep("  " .. frag, D-1)
end

-- Solve the problem
solve(mat, sol, foundsol)
