export rho_init,
        mpdo_I,
        measure_mpdo_single,
        measure_mpdo_long


#--------------Define Constants----------------#
const sz = (1/2.0)*[1. 0.;0. -1.]
const sx = (1/2.0)*[0. 1.; 1. 0.]
const sy = (1/2.0)*[0. -1.0im; +1.0im 0.]
const sp = [0. 1. ; 0. 0.]
const sm = [0. 0. ; 1. 0.]
const id_half = [1. 0.; 0. 1.]
const Pdn = [0. 0. ; 0. 1.]
const Pup = [1. 0. ; 0. 0.]


const sz_1 = [1. 0. 0.;0. 0. 0.;0. 0. -1.]
const sx_1 = (1.0/sqrt(2.0))*[0. 1. 0.;1. 0. 1.;0. 1. 0.]
const sy_1 = (-1.0im/sqrt(2.0))*[0. 1. 0.;-1. 0. 1.;0. -1. 0.]
const sp_1 = sqrt(2.0)*[0. 1. 0;0. 0. 1.;0. 0. 0.]
const sm_1 = sqrt(2.0)*[0. 0. 0;1. 0. 0.;0. 1. 0.]
const id_1 = [1. 0. 0.; 0. 1. 0.; 0. 0. 1.]

# const pauli = Dict("I"=>Matrix(1.0I,2,2) ,"sx"=> sx,
#                    "sy"=>sy, "sz"=>sz, "sp"=>sp, "sm"=>sm,
#                    "sPup"=> Pup, "sPdn"=>Pdn)


#-------------Define Site types-------------------#

ITensors.space(::SiteType"mpdo=1") = 9
ITensors.space(::SiteType"mpdo=1/2") = 4

#--------------Define Ops for spin=1/2-------------------#




function ITensors.op!(Op::ITensor,O::OpName"TPdn",::SiteType"mpdo=1/2",s::Index)
    A=kron(id_half, Pdn)
    for i in 1:4
        for j in 1:4
            Op[s'=>i,s=>j] = A[i,j]
        end
    end
end
        
        
function ITensors.op!(Op::ITensor,O::OpName"SPdn",::SiteType"mpdo=1/2",s::Index)
    A=kron(Pdn, id_half)
    for i in 1:4
        for j in 1:4
            Op[s'=>i,s=>j] = A[i,j]
        end
    end
end
        
        
function ITensors.op!(Op::ITensor,O::OpName"SPup",::SiteType"mpdo=1/2",s::Index)
    A=kron(Pup, id_half)
    for i in 1:4
        for j in 1:4
            Op[s'=>i,s=>j] = A[i,j]
        end
    end
end
        
        
function ITensors.op!(Op::ITensor,O::OpName"TPup",::SiteType"mpdo=1/2",s::Index)
    A=kron(id_half, Pup)
    for i in 1:4
        for j in 1:4
            Op[s'=>i,s=>j] = A[i,j]
        end
    end
end

function ITensors.op!(Op::ITensor,O::OpName"Tz",::SiteType"mpdo=1/2",s::Index)
    A=kron(id_half, sz)
    for i in 1:4
        for j in 1:4
            Op[s'=>i,s=>j] = A[i,j]
        end
    end
end

function ITensors.op!(Op::ITensor,O::OpName"Tx",::SiteType"mpdo=1/2",s::Index)
    A=kron(id_half, sx)
    for i in 1:4
        for j in 1:4
            Op[s'=>i,s=>j] = A[i,j]
        end
    end
end

function ITensors.op!(Op::ITensor,O::OpName"Ty",::SiteType"mpdo=1/2",s::Index)
    A=kron(id_half, sy)
    for i in 1:4
        for j in 1:4
            Op[s'=>i,s=>j] = A[i,j]
        end
    end
end


function ITensors.op!(Op::ITensor,O::OpName"Sz",::SiteType"mpdo=1/2",s::Index)
    A=kron(sz, id_half)
    for i in 1:4
        for j in 1:4
            Op[s'=>i,s=>j] = A[i,j]
        end
    end
end

function ITensors.op!(Op::ITensor,O::OpName"Sx",::SiteType"mpdo=1/2",s::Index)
    A=kron(sx, id_half)
    for i in 1:4
        for j in 1:4
            Op[s'=>i,s=>j] = A[i,j]
        end
    end
end

function ITensors.op!(Op::ITensor,O::OpName"Sy",::SiteType"mpdo=1/2",s::Index)
    A=kron(sy, id_half)
    for i in 1:4
        for j in 1:4
            Op[s'=>i,s=>j] = A[i,j]
        end
    end
end

function ITensors.op!(Op::ITensor,O::OpName"SpTp",::SiteType"mpdo=1/2",s::Index)
    A=kron(sp, sp)
    for i in 1:4
        for j in 1:4
            Op[s'=>i,s=>j] = A[i,j]
        end
    end
end


function ITensors.op!(Op::ITensor,O::OpName"SmTm",::SiteType"mpdo=1/2",s::Index)
    A=kron(sm, sm)
    for i in 1:4
        for j in 1:4
            Op[s'=>i,s=>j] = A[i,j]
        end
    end
end

function ITensors.op!(Op::ITensor,O::OpName"SpSm",::SiteType"mpdo=1/2",s::Index)
    A=kron(*(sp,sm), id_half)
    for i in 1:4
        for j in 1:4
            Op[s'=>i,s=>j] = A[i,j]
        end
    end
end


function ITensors.op!(Op::ITensor,O::OpName"SmSp",::SiteType"mpdo=1/2",s::Index)
    A=kron(*(sm,sp), id_half)
    for i in 1:4
        for j in 1:4
            Op[s'=>i,s=>j] = A[i,j]
        end
    end
end


function ITensors.op!(Op::ITensor,O::OpName"TmTp",::SiteType"mpdo=1/2",s::Index)
    A=kron(id_half, *(sm,sp))
    for i in 1:4
        for j in 1:4
            Op[s'=>i,s=>j] = A[i,j]
        end
    end
end


function ITensors.op!(Op::ITensor,O::OpName"TpTm",::SiteType"mpdo=1/2",s::Index)
    A=kron(id_half, *(sp,sm))
    for i in 1:4
        for j in 1:4
            Op[s'=>i,s=>j] = A[i,j]
        end
    end
end



function ITensors.op!(Op::ITensor,O::OpName"Tp",::SiteType"mpdo=1/2",s::Index)
    A=kron(id_half, sp)
    for i in 1:4
        for j in 1:4
            Op[s'=>i,s=>j] = A[i,j]
        end
    end
end
        
function ITensors.op!(Op::ITensor,O::OpName"Tm",::SiteType"mpdo=1/2",s::Index)
    A=kron(id_half, sm)
    for i in 1:4
        for j in 1:4
            Op[s'=>i,s=>j] = A[i,j]
        end
    end
end
        
        
function ITensors.op!(Op::ITensor,O::OpName"Sm",::SiteType"mpdo=1/2",s::Index)
    A=kron(sm, id_half)
    for i in 1:4
        for j in 1:4
            Op[s'=>i,s=>j] = A[i,j]
        end
    end
end
        
function ITensors.op!(Op::ITensor,O::OpName"Sp",::SiteType"mpdo=1/2",s::Index)
    A=kron(sp, id_half)
    for i in 1:4
        for j in 1:4
            Op[s'=>i,s=>j] = A[i,j]
        end
    end
end


#-------------Define Ops for Spin=1-----------------#

function ITensors.op!(Op::ITensor,O::OpName"Tz",::SiteType"mpdo=1",s::Index)
    A=kron(id_1, sz_1)
    for i in 1:9
        for j in 1:9
            Op[s'=>i,s=>j] = A[i,j]
        end
    end
end

function ITensors.op!(Op::ITensor,O::OpName"Tx",::SiteType"mpdo=1",s::Index)
    A=kron(id_1, sx_1)
    for i in 1:9
        for j in 1:9
            Op[s'=>i,s=>j] = A[i,j]
        end
    end
end

function ITensors.op!(Op::ITensor,O::OpName"Ty",::SiteType"mpdo=1",s::Index)
    A=kron(id_1, sy_1)
    for i in 1:9
        for j in 1:9
            Op[s'=>i,s=>j] = A[i,j]
        end
    end
end


function ITensors.op!(Op::ITensor,O::OpName"Sz",::SiteType"mpdo=1",s::Index)
    A=kron(sz_1, id_1)
    for i in 1:9
        for j in 1:9
            Op[s'=>i,s=>j] = A[i,j]
        end
    end
end

function ITensors.op!(Op::ITensor,O::OpName"Sx",::SiteType"mpdo=1",s::Index)
    A=kron(sx_1, id_1)
    for i in 1:9
        for j in 1:9
            Op[s'=>i,s=>j] = A[i,j]
        end
    end
end

function ITensors.op!(Op::ITensor,O::OpName"Sy",::SiteType"mpdo=1",s::Index)
    A=kron(sy_1, id_1)
    for i in 1:9
        for j in 1:9
            Op[s'=>i,s=>j] = A[i,j]
        end
    end
end

function ITensors.op!(Op::ITensor,O::OpName"SpTp",::SiteType"mpdo=1",s::Index)
    A=kron(sp_1, sp_1)
    for i in 1:9
        for j in 1:9
            Op[s'=>i,s=>j] = A[i,j]
        end
    end
end


function ITensors.op!(Op::ITensor,O::OpName"SmTm",::SiteType"mpdo=1",s::Index)
    A=kron(sm_1, sm_1)
    for i in 1:9
        for j in 1:9
            Op[s'=>i,s=>j] = A[i,j]
        end
    end
end

function ITensors.op!(Op::ITensor,O::OpName"SpSm",::SiteType"mpdo=1",s::Index)
    A=kron(*(sp_1,sm_1), id_1)
    for i in 1:9
        for j in 1:9
            Op[s'=>i,s=>j] = A[i,j]
        end
    end
end


function ITensors.op!(Op::ITensor,O::OpName"SmSp",::SiteType"mpdo=1",s::Index)
    A=kron(*(sm_1,sp_1), id_1)
    for i in 1:9
        for j in 1:9
            Op[s'=>i,s=>j] = A[i,j]
        end
    end
end


function ITensors.op!(Op::ITensor,O::OpName"TmTp",::SiteType"mpdo=1",s::Index)
    A=kron(id_1, *(sm_1,sp_1))
    for i in 1:9
        for j in 1:9
            Op[s'=>i,s=>j] = A[i,j]
        end
    end
end


function ITensors.op!(Op::ITensor,O::OpName"TpTm",::SiteType"mpdo=1",s::Index)
    A=kron(id_1, *(sp_1,sm_1))
    for i in 1:9
        for j in 1:9
            Op[s'=>i,s=>j] = A[i,j]
        end
    end
end

function ITensors.op!(Op::ITensor,O::OpName"Sxx",::SiteType"mpdo=1",s::Index)
    A=kron(*(sx_1, sx_1), id_1)
    for i in 1:9
        for j in 1:9
            Op[s'=>i,s=>j] = A[i,j]
        end
    end
end
        
function ITensors.op!(Op::ITensor,O::OpName"Sxy",::SiteType"mpdo=1",s::Index)
    A=kron(*(sx_1, sy_1), id_1)
    for i in 1:9
        for j in 1:9
            Op[s'=>i,s=>j] = A[i,j]
        end
    end
end
        
        
function ITensors.op!(Op::ITensor,O::OpName"Sxz",::SiteType"mpdo=1",s::Index)
    A=kron(*(sx_1, sz_1), id_1)
    for i in 1:9
        for j in 1:9
            Op[s'=>i,s=>j] = A[i,j]
        end
    end
end
     
function ITensors.op!(Op::ITensor,O::OpName"Syx",::SiteType"mpdo=1",s::Index)
    A=kron(*(sy_1, sx_1), id_1)
    for i in 1:9
        for j in 1:9
            Op[s'=>i,s=>j] = A[i,j]
        end
    end
end
        
        
function ITensors.op!(Op::ITensor,O::OpName"Syy",::SiteType"mpdo=1",s::Index)
    A=kron(*(sy_1, sy_1), id_1)
    for i in 1:9
        for j in 1:9
            Op[s'=>i,s=>j] = A[i,j]
        end
    end
end
        
function ITensors.op!(Op::ITensor,O::OpName"Syz",::SiteType"mpdo=1",s::Index)
    A=kron(*(sy_1, sz_1), id_1)
    for i in 1:9
        for j in 1:9
            Op[s'=>i,s=>j] = A[i,j]
        end
    end
end
        
        
function ITensors.op!(Op::ITensor,O::OpName"Szx",::SiteType"mpdo=1",s::Index)
    A=kron(*(sz_1, sx_1), id_1)
    for i in 1:9
        for j in 1:9
            Op[s'=>i,s=>j] = A[i,j]
        end
    end
end
        
function ITensors.op!(Op::ITensor,O::OpName"Szy",::SiteType"mpdo=1",s::Index)
    A=kron(*(sz_1, sy_1), id_1)
    for i in 1:9
        for j in 1:9
            Op[s'=>i,s=>j] = A[i,j]
        end
    end
end
        
function ITensors.op!(Op::ITensor,O::OpName"Szz",::SiteType"mpdo=1",s::Index)
    A=kron(*(sz_1, sz_1), id_1)
    for i in 1:9
        for j in 1:9
            Op[s'=>i,s=>j] = A[i,j]
        end
    end
end
        
function ITensors.op!(Op::ITensor,O::OpName"Txx",::SiteType"mpdo=1",s::Index)
    A=kron(id_1, *(sx_1, sx_1))
    for i in 1:9
        for j in 1:9
            Op[s'=>i,s=>j] = A[i,j]
        end
    end
end
        
function ITensors.op!(Op::ITensor,O::OpName"Txy",::SiteType"mpdo=1",s::Index)
    A=kron(id_1, *(sx_1, sy_1))
    for i in 1:9
        for j in 1:9
            Op[s'=>i,s=>j] = A[i,j]
        end
    end
end
        
function ITensors.op!(Op::ITensor,O::OpName"Txz",::SiteType"mpdo=1",s::Index)
    A=kron(id_1, *(sx_1, sz_1))
    for i in 1:9
        for j in 1:9
            Op[s'=>i,s=>j] = A[i,j]
        end
    end
end
        
        
function ITensors.op!(Op::ITensor,O::OpName"Tyx",::SiteType"mpdo=1",s::Index)
    A=kron(id_1, *(sy_1, sx_1))
    for i in 1:9
        for j in 1:9
            Op[s'=>i,s=>j] = A[i,j]
        end
    end
end
        
function ITensors.op!(Op::ITensor,O::OpName"Tyy",::SiteType"mpdo=1",s::Index)
    A=kron(id_1, *(sy_1, sy_1))
    for i in 1:9
        for j in 1:9
            Op[s'=>i,s=>j] = A[i,j]
        end
    end
end
        
function ITensors.op!(Op::ITensor,O::OpName"Tyz",::SiteType"mpdo=1",s::Index)
    A=kron(id_1, *(sy_1, sz_1))
    for i in 1:9
        for j in 1:9
            Op[s'=>i,s=>j] = A[i,j]
        end
    end
end
        
function ITensors.op!(Op::ITensor,O::OpName"Tzx",::SiteType"mpdo=1",s::Index)
    A=kron(id_1, *(sz_1, sx_1))
    for i in 1:9
        for j in 1:9
            Op[s'=>i,s=>j] = A[i,j]
        end
    end
end
        
function ITensors.op!(Op::ITensor,O::OpName"Tzy",::SiteType"mpdo=1",s::Index)
    A=kron(id_1, *(sz_1, sy_1))
    for i in 1:9
        for j in 1:9
            Op[s'=>i,s=>j] = A[i,j]
        end
    end
end
        
function ITensors.op!(Op::ITensor,O::OpName"Tzz",::SiteType"mpdo=1",s::Index)
    A=kron(id_1, *(sz_1, sz_1))
    for i in 1:9
        for j in 1:9
            Op[s'=>i,s=>j] = A[i,j]
        end
    end
end




#-----------Define basic density matrices-------------------#

function Up_site_rho(d::Int)
    if (d==4)
        return [1. 0.;0. 0.]
    end
    if (d==9)
        return [1. 0. 0.;0. 0. 0.;0. 0. 0.]
    else 
        throw("correct your site type")    
    end
end


function Dn_site_rho(d::Int)
    
    if (d==4)
        return [0. 0.;0. 1.]
    end
    if (d==9)
        return [0. 0. 0.;0. 0. 0.;0. 0. 1.]
    else 
        throw("correct your site type")    
    end
end

function Identity_site_rho(d::Int)
    if (d==4)
        mx = 1.0/sqrt(2.)
        return mx*[1. 0.;0. 1.]
    end
    if (d==9)
        mx = 1.0/sqrt(3.)
        return mx*[1. 0. 0.;0. 1. 0.;0. 0. 1.]
     else 
        throw("correct your site type")    
    end
end

# function rand_single_spin_wf()
#     c = rand(ComplexF64,2)
#     return c/LinearAlgebra.norm(c)
# end


#---------general definition of which_T function-------------------#

function which_T(which,j::Int, d::Int)
    if which==:allup
        return Up_site_rho(d)
    elseif which==:alldown
        return Dn_site_rho(d)
    elseif which==:I  
        return Identity_site_rho(d)
    elseif which==:afm 
        return j%2==1 ? Up_site_rho(d) : Dn_site_rho(d)    
    #elseif which==:rand_pure
    #    c = rand_single_spin_wf(); 
    #    return c*c'
    else 
        throw("initial state $which is invalid")    
    end
end


#---------------Initialize rho--------------------------#

function rho_init(sites::Vector{<:Index},linkdim::Int,which=:allup;kwargs...)   
    N = length(sites)
    M =  MPS(sites)
    d = dim(sites[N])

    if N==1
        M[1] = ITensor(reshape(which_T(which,1,d), (d)), sites[1]) 
        M[1] /= norm(M[1])
        return M
    end

    l = Vector{Index}(undef,N) 
    chi = linkdim 
    
    l[N-1] = Index(chi,"Link,l=$(N-1)")
    T = reshape(which_T(which,N,d), (chi,d))
    M[N] = ITensor(T, l[N-1],sites[N])

    for j=N-1:-1:2
        l[j-1] = Index(chi,"Link,l=$(j-1)")
        T = reshape(which_T(which,j,d), (chi,d,chi))
        M[j] = ITensor(T, l[j-1],sites[j],l[j])
    end
    
    l0 = Index(1,"Link,l=0")
    T = reshape(which_T(which,1,d), (d,chi))
    M[1] = ITensor(T, l0,sites[1],l[1])
    M[1] *= setelt(l0=>1)

    M.llim = 0
    M.rlim = 2  #what are these two?
    
    return M
end


#-----------Rho for infinite temp-----------------------#

function mpdo_I(psi::MPS)
    s = siteinds(psi)
    N = length(s)
    M =  MPS(s)
    d = dim(s[1])
    
    if (d==4)
        matID = id_half
    else (d==9)
        matID = id_1
    end
    
    #l = linkinds(psi) 
    l = Vector{Index}(undef,N) 
    chi = 1 
    
    l[N-1] = Index(chi,"Link,l=$(N-1)")
    T = reshape(matID, (chi,d))
    M[N] = ITensor(T, l[N-1],s[N])

    for j=N-1:-1:2
        l[j-1] = Index(chi,"Link,l=$(j-1)")
        T = reshape(matID, (chi,d,chi))
        M[j] = ITensor(T, l[j-1],s[j],l[j])
    end
    
    l0 = Index(1,"Link,l=0")
    T = reshape(matID, (d,chi))
    M[1] = ITensor(T, l0,s[1],l[1])
    M[1] *= setelt(l0=>1)

    M.llim = 0
    M.rlim = 2  #what are these two?
    
    return M
end

function measure_mpdo_single(psi::MPS, O::String, i::Int)   #mdpo_I as input
    phi=deepcopy(psi)
    phi[i]=noprime(op(O,siteind(psi,i))*phi[i])
    return inner(mpdo_I(phi),phi)/inner(mpdo_I(psi),psi)
end

function measure_mpdo_long(psi::MPS, Os::Vector{String}, is::Vector{Int})
    phi=deepcopy(psi)
    
    if size(is,1)==1
        return measure_mpdo_single(psi,Os[1],is[1])
    end

    for (k,O) in enumerate(Os)
        phi[is[k]]=noprime(op(O,siteind(psi,is[k]))*phi[is[k]])
    end

    return inner(mpdo_I(phi),phi)/inner(mpdo_I(psi),psi)
end

# function measure_long(psi::MPS, Os::Vector{String}, is::Vector{Int})
#     phi=deepcopy(psi)
    
#     if size(is,1)==1
#         return measure_mpdo_single(psi,Os[1],is[1])
#     end

#     for (k,O) in enumerate(Os)
#         phi[is[k]]=noprime(op(O,siteind(psi,is[k]))*phi[is[k]])
#     end

#     return inner(mpdo_I(phi),phi)/inner(mpdo_I(psi),psi)
# end


# function measure_single(psi::MPS, O::String, i::Int)	#mdpo_I as input
#     phi=deepcopy(psi)
#     phi[i]=noprime(op(O,siteind(psi,i))*phi[i])
#     return inner(mpdo_I(phi),phi)/inner(mpdo_I(psi),psi)
# end
