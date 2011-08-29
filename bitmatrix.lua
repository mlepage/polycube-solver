-- Bit matrix module
-- Copyright (C) 2009 Marc A. Lepage

require "bit"
local band, bor, bnot = bit.band, bit.bor, bit.bnot
local lshift, rshift = bit.lshift, bit.rshift
local assert, ceil, write = assert, math.ceil, io.write

module(...)

-- Matrix bit group size
local G = 32

-- Masks for bit operations
local mask, mask_n, mask_lt, mask_ge, mask_gt, mask_o, mask_u = {}, {}, {}, {}, {}, lshift(1, G-1), 1
for n = 0, G-1 do
	mask[n]    = lshift(1, n)                -- only bit n is set
	mask_n[n]  = bnot(mask[n])               -- all bits but n are set
	mask_lt[n] = mask[n] - 1                 -- all bits less than n are set
	mask_ge[n] = bnot(mask_lt[n])            -- all bits greater than or equal to n are set
	mask_gt[n] = band(mask_ge[n], mask_n[n]) -- all bits greater than n are set
--	mask_o                                   -- only bit G is set (overflow)
--	mask_u                                   -- only bit 1 is set (underflow)
end

-- Creates an MxN matrix with M rows and N columns
function create(m, n)
	m, n = m and m or 0, n and n or 0
	assert(0 <= m and 0 <= n)
	local mat = { m=m, n=n }
	for i = 1, mat.m do
		mat[i] = {}
		for k = 1, ceil(mat.n/G) do
			mat[i][k] = 0
		end
	end
	return mat
end

-- Copies a matrix
function copy(mat)
	local matc = { m=mat.m, n=mat.n}
	for i = 1, mat.m do
		matc[i] = {}
		for k = 1, ceil(mat.n/G) do
			matc[i][k] = mat[i][k]
		end
	end
	return matc
end

function get(mat, i, j)
	assert(1 <= i and i <= mat.m and 1 <= j and j <= mat.n)
	local k, l = ceil(j/G), (j-1)%G
	return band(mat[i][k], mask[l]) ~= 0 and 1 or 0
end

function set(mat, i, j, v)
	assert(1 <= i and i <= mat.m and 1 <= j and j <= mat.n)
	local k, l = ceil(j/G), (j-1)%G
	if v ~= 0 then
		mat[i][k] = bor(mat[i][k], mask[l])
	else
		mat[i][k] = band(mat[i][k], mask_n[l])
	end
end

function insertrow(mat, i)
	assert(1 <= i and i <= mat.m+1)
	for i = mat.m, i, -1 do
		mat[i+1] = mat[i]
	end
	mat[i] = {}
	for k = 1, ceil(mat.n/G) do
		mat[i][k] = 0
	end
	mat.m = mat.m+1
end

function removerow(mat, i)
	assert(1 <= i and i <= mat.m)
	for i = i, mat.m-1 do
		mat[i] = mat[i+1]
	end
	mat[mat.m] = nil
	mat.m = mat.m-1
end

function insertcol(mat, j)
	assert(1 <= j and j <= mat.n+1)
	local kn, k, l = ceil(mat.n/G), ceil(j/G), (j-1)%G
	local mask_lt, mask_ge = mask_lt[l], mask_ge[l]
	local grow = kn < ceil((mat.n+1)/G)
	for i = 1, mat.m do
		if grow then mat[i][kn+1] = 0 end
		local o = band(mask_o, mat[i][k]) ~= 0 and mask_u or 0
		mat[i][k] = bor(lshift(band(mask_ge, mat[i][k]), 1), band(mask_lt, mat[i][k]))
		for k = k+1, kn do
			o, mat[i][k] = band(mask_o, mat[i][k]) ~= 0 and mask_u or 0, bor(lshift(mat[i][k], 1), o)
		end
	end
	mat.n = mat.n+1
end

function removecol(mat, j)
	assert(1 <= j and j <= mat.n)
	local kn, k, l = ceil(mat.n/G), ceil(j/G), (j-1)%G
	local mask_lt, mask_gt = mask_lt[l], mask_gt[l]
	local shrink = ceil((mat.n-1)/G) < kn
	for i = 1, mat.m do
		local u = 0
		for k = kn, k+1, -1 do
			u, mat[i][k] = band(mask_u, mat[i][k]) ~= 0 and mask_o or 0, bor(u, rshift(mat[i][k], 1))
		end
		mat[i][k] = bor(u, rshift(band(mask_gt, mat[i][k]), 1), band(mask_lt, mat[i][k]))
		if shrink then mat[i][kn] = nil end
	end
	mat.n = mat.n-1
end

function print(mat)
	for i = 1, mat.m do
		for j = 1, mat.n do
			write(get(mat, i, j))
		end
		write("\n")
	end
end
