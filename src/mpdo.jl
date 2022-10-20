export rho_init,
        mpdo_I,
        measure_mpdo_single,
        measure_mpdo_long

# const sz = [1. 0.;0. -1.]
# const sx = [0. 1.; 1. 0.]
# const sy = [0. -1.0im; +1.0im 0.]
# const sp = [0. 1. ; 0. 0.]
# const sm = [0. 0. ; 1. 0.]
# const Pup = 1/2*(I(2)+sz)
# const Pdn = 1/2*(I(2)-sz)

# const pauli = Dict("I"=>Matrix(1.0I,2,2) ,"sx"=> sx,
#                    "sy"=>sy, "sz"=>sz, "sp"=>sp, "sm"=>sm,
#                    "sPup"=> Pup, "sPdn"=>Pdn)

ITensors.space(::SiteType"mpdo") = 4

function ITensors.op!(Op::ITensor,		
                      O::OpName"Tx",
                      ::SiteType"mpdo",
                      s::Index)

	Op[s'=>1,s=>2]=1
	Op[s'=>2,s=>1]=1
	Op[s'=>3,s=>4]=1
	Op[s'=>4,s=>3]=1
    # A=kron(I(2), sx)
    # for i in 1:4
    #     for j in 1:4
    #         Op[s'=>i,s=>j] = A[i,j]
    #     end
    # end
end

function ITensors.op!(Op::ITensor,
                      O::OpName"Ty",
                      ::SiteType"mpdo",
                      s::Index)
	ITensors.complex!(Op)
	Op[s'=>1,s=>2]=-1.0im
	Op[s'=>2,s=>1]=1.0im
	Op[s'=>3,s=>4]=-1.0im
	Op[s'=>4,s=>3]=1.0im
end

function ITensors.op!(Op::ITensor,
                      O::OpName"Tz",
                      ::SiteType"mpdo",
                      s::Index)
    Op[s'=>1,s=>1]=1
    Op[s'=>2,s=>2]=-1
    Op[s'=>3,s=>3]=1
    Op[s'=>4,s=>4]=-1
end

function ITensors.op!(Op::ITensor,
                      O::OpName"Id",
                      ::SiteType"mpdo",
                      s::Index)
    Op[s'=>1,s=>1]=1
    Op[s'=>2,s=>2]=1
    Op[s'=>3,s=>3]=1
    Op[s'=>4,s=>4]=1
end

function ITensors.op!(Op::ITensor,
                      O::OpName"Sx",
                      ::SiteType"mpdo",
                      s::Index)
    Op[s'=>1,s=>3]=1
    Op[s'=>2,s=>4]=1
    Op[s'=>3,s=>1]=1
    Op[s'=>4,s=>2]=1
end

function ITensors.op!(Op::ITensor,
                      O::OpName"Sy",
                      ::SiteType"mpdo",
                      s::Index)
                      #::Type{T})
    ITensors.complex!(Op)
    Op[s'=>1,s=>3]=-1.0im
    Op[s'=>2,s=>4]=-1.0im
    Op[s'=>3,s=>1]=1.0im
    Op[s'=>4,s=>2]=1.0im
end

function ITensors.op!(Op::ITensor,
                      O::OpName"Sz",
                      ::SiteType"mpdo",
                      s::Index)
    Op[s'=>1,s=>1]=1
    Op[s'=>2,s=>2]=1
    Op[s'=>3,s=>3]=-1
    Op[s'=>4,s=>4]=-1
end

function ITensors.op!(Op::ITensor,
                      O::OpName"Sp",
                      ::SiteType"mpdo",
                      s::Index)
    Op[s'=>1,s=>3]=1
    Op[s'=>2,s=>4]=1
end

function ITensors.op!(Op::ITensor,
                      O::OpName"Sm",
                      ::SiteType"mpdo",
                      s::Index)
    Op[s'=>3,s=>1]=1
    Op[s'=>4,s=>2]=1
end

function ITensors.op!(Op::ITensor,
                      O::OpName"Tp",
                      ::SiteType"mpdo",
                      s::Index)
    Op[s'=>1,s=>2]=1
    Op[s'=>3,s=>4]=1
end

function ITensors.op!(Op::ITensor,
                      O::OpName"Tm",
                      ::SiteType"mpdo",
                      s::Index)
    Op[s'=>2,s=>1]=1
    Op[s'=>4,s=>3]=1
end

function ITensors.op!(Op::ITensor,      
                      O::OpName"SpTp",
                      ::SiteType"mpdo",
                      s::Index)
    # A=kron(sp, sp)
    # for i in 1:4
    #     for j in 1:4
    #         Op[s'=>i,s=>j] = A[i,j]
    #     end
    # end
    Op[s'=>1,s=>4]=1
end

function ITensors.op!(Op::ITensor,      
                      O::OpName"SmTm",
                      ::SiteType"mpdo",
                      s::Index)
    # A=kron(sm, sm)
    # for i in 1:4
    #     for j in 1:4
    #         Op[s'=>i,s=>j] = A[i,j]
    #     end
    # end
    Op[s'=>4,s=>1]=1
end

function ITensors.op!(Op::ITensor,      
                      O::OpName"SmSp",
                      ::SiteType"mpdo",
                      s::Index)
    # A=kron(Pdn, I(2))
    # for i in 1:4
    #     for j in 1:4
    #         Op[s'=>i,s=>j] = A[i,j]
    #     end
    # end
    Op[s'=>3,s=>3]=1
    Op[s'=>4,s=>4]=1
end

function ITensors.op!(Op::ITensor,      
                      O::OpName"SpSm",
                      ::SiteType"mpdo",
                      s::Index)
    # A=kron(Pup, I(2))
    # for i in 1:4
    #     for j in 1:4
    #         Op[s'=>i,s=>j] = A[i,j]
    #     end
    # end
    Op[s'=>1,s=>1]=1
    Op[s'=>2,s=>2]=1
end

function ITensors.op!(Op::ITensor,      
                      O::OpName"TmTp",
                      ::SiteType"mpdo",
                      s::Index)
    # A=kron(I(2),Pdn)
    # for i in 1:4
    #     for j in 1:4
    #         Op[s'=>i,s=>j] = A[i,j]
    #     end
    # end
    Op[s'=>2,s=>2]=1
    Op[s'=>4,s=>4]=1
end

function ITensors.op!(Op::ITensor,      
                      O::OpName"TpTm",
                      ::SiteType"mpdo",
                      s::Index)
    # A=kron(I(2),Pup)
    # for i in 1:4
    #     for j in 1:4
    #         Op[s'=>i,s=>j] = A[i,j]
    #     end
    # end
    Op[s'=>1,s=>1]=1
    Op[s'=>3,s=>3]=1
end

const up_site_rho = [1.0+0.0im 0; 0. 0.]    #ZL make them function of sitetype
const dn_site_rho = [0. 0.; 0. 1.0+0im]
const Id_site_rho = Matrix((1.0+0.0im)*I,2,2)/sqrt(2)

function rand_single_spin_wf()
    c = rand(ComplexF64,2)
    return c/LinearAlgebra.norm(c)
end

function which_T(which,j::Int)
    if which==:allup
        return up_site_rho
    elseif which==:alldown
        return dn_site_rho
    elseif which==:I  
        return Id_site_rho
    elseif which==:afm 
        return j%2==1 ? up_site_rho : dn_site_rho    
    elseif which==:rand_pure
        c = rand_single_spin_wf(); 
        return c*c'
    else 
        throw("initial state $which is invalid")    
    end
end

function rho_init(sites::Vector{<:Index},linkdim::Int,which=:allup;kwargs...)   
    N = length(sites)
    M =  MPS(sites)
    d = dim(sites[N])

    if N==1
        M[1] = ITensor(reshape(which_T(which,1), (d)), sites[1]) 
        M[1] /= norm(M[1])
        return M
    end

    l = Vector{Index}(undef,N) 
    chi = linkdim 
    
    l[N-1] = Index(chi,"Link,l=$(N-1)")
    T = reshape(which_T(which,N), (chi,d))
    M[N] = ITensor(T, l[N-1],sites[N])

    for j=N-1:-1:2
        l[j-1] = Index(chi,"Link,l=$(j-1)")
        T = reshape(which_T(which,j), (chi,d,chi))
        M[j] = ITensor(T, l[j-1],sites[j],l[j])
    end
    
    l0 = Index(1,"Link,l=0")
    T = reshape(which_T(which,1), (d,chi))
    M[1] = ITensor(T, l0,sites[1],l[1])
    M[1] *= setelt(l0=>1)

    M.llim = 0
    M.rlim = 2  #what are these two?
    
    return M
end

function mpdo_I(psi::MPS)
    s = siteinds(psi)
    N = length(s)
    M =  MPS(s)
    d = dim(s[1])
    
    matID=Matrix((1.0+0.0im)*I,2,2)
    
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

function measure_mpdo_single(psi::MPS, O::String, i::Int)	#mdpo_I as input
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


# ---------------------------------------------
# I should check this code again

"""
function reduced_rho(psi::MPS, i1::Int, i2::Int)
reduced density matrix between sites i1 and i2
psi is the input density matrix
"""
function reduced_rho(psi::MPS, i1::Int, i2::Int)
    
    if (i2-i1+1 > length(siteinds(psi))) || i1>i2
        throw("invalid range") 
    end
    
    #phi=deepcopy(psi) #why not needed if I don;'t do any reorth'
    rho_Id=mpdo_I(psi)
    N=length(siteinds(psi))
    
    #left environment
    L=psi[1]*rho_Id[1]
    for j in 2:i1-1
        L*=psi[j]*rho_Id[j]
    end
    
    #right environment
    R=psi[end]*rho_Id[end]
    for j in N-1:-1:i2+1
        R*=psi[j]*rho_Id[j]
    end
    
    M = MPS(i2-i1+1)
    for j in i1:i2
        M[j-i1+1] = psi[j]
    end
    
    M[i2-i1+1]*=R
    M[1]*=L
    
    return M
end



"""
opeentropy(psi::MPS, b::Int)
calculated operator entanglement entropy for rho=psi
at site b
"""
function opeentropy(psi::MPS, b::Int)
    orthogonalize!(psi, b)  #why do I need ortogonalize?
    U,S,V = svd(psi[b], (linkind(psi, b-1), siteind(psi,b)))
    SvN = 0.0
    norma= 0.0
    for n=1:dim(S, 1)
        p = S[n,n]^2
        SvN -= p * log(p)
        norma+=p
    end
    norma≈1.0 ? nothing : println("not normalized, sum(S_i^2)=",norma)
    return SvN
end

"""
checkdone!(cb,i,deriv_tol)
calculates running average of the derivative of the average of operators 
measured in cb::LocalMeasurementCallback at site i and returns true if 
the running derivate is smaller than deriv_tol
"""
function checkdone!(cb,i,deriv_tol)
    
    datT=getindex.(measurements(cb)[cb.ops[1]],i)
    for ik in 2:length(cb.ops)
        datT = datT .+ getindex.(measurements(cb)[cb.ops[ik]],i)
    end
    datT./length(cb.ops)
    
    if length(cb.ts)<10
        return false
    else
        dM = [(datT[i+1] - datT[i])/dt
                for i in length(datT)-10:length(datT)-1]
        dM = mean(abs.(dM))
        dM<deriv_tol ? (return true) : (return false)
    end
end