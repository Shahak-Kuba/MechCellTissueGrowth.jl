Xᵩ(T) = T
Yᵩ(T) = 2 + 0.5*cos(3*T)

# Circular Boundary
X(R,θ) = R.*cos.(θ);
Y(R,θ) = R.*sin.(θ);

# for all other regular regular polygons

function polygon_vertices(N, R, theta)
    # N: Number of sides
    # R: Radius
    # theta: Initial rotation angle in radians

    vertices = []
    for i in 0:(N-1)
        angle = theta + (2 * π * i / N)
        x = R * cos(angle)
        y = R * sin(angle)
        push!(vertices, (x, y))
    end
    return vertices
end

function interpolate_segment(p1, p2, w, dist_type)
    [((1 - t) .* p1 .+ t .* p2) for t in NodeDistribution(0, 1, w, dist_type)]
end

# for irregular polygons
function CrossVertecies(side_length, offset)
    CrossVerts = []
    push!(CrossVerts,(side_length + offset, offset))
    push!(CrossVerts,(offset, offset))
    push!(CrossVerts,(offset, offset + side_length))
    push!(CrossVerts,(-offset, offset + side_length))
    push!(CrossVerts,(-offset, offset))
    push!(CrossVerts,(-offset - side_length, offset))
    push!(CrossVerts,(-offset - side_length, -offset))
    push!(CrossVerts,(-offset, -offset))
    push!(CrossVerts,(-offset, -offset - side_length))
    push!(CrossVerts,(offset, -offset - side_length))
    push!(CrossVerts,(offset, -offset))
    push!(CrossVerts,(offset + side_length, -offset))
    #push!(CrossVerts,(side_length + offset, offset))

    #return hcat(CrossVerts...)'
    return CrossVerts
end

function StarVerticies(N, R, Rotation_Angle, rotation_angle)
    # finding R such that areas will match with formula A = 2N(0.5*R₀*rₒ*sin(θ))
    R₀ = √((4*π*(R^2))/(2*N*sin(π/N)))
    rₒ = R₀/2
    # empty vector
    StarVerts = []
    # generating verticies for outside polygon
    VERTS = regular_polygon_vertices(N, R₀, Rotation_Angle)
    # generating verticies for inside polygon
    verts = regular_polygon_vertices(N, rₒ, rotation_angle)

    # combining the two
    for i in 1:N
        push!(StarVerts, (VERTS[i,1], VERTS[i,2]))
        push!(StarVerts, (verts[i,1], verts[i,2]))
    end
    #push!(StarVerts, VERTS[1,:])

    #return hcat(StarVerts...)'
    return StarVerts
end

function regular_polygon_vertices(N, R, rotation_angle)
    vertices = Vector{Float64}[]

    for i in 0:N-1
        angle = 2π * i / N + rotation_angle
        x = R * cos(angle)
        y = R * sin(angle)
        push!(vertices, [x, y])
    end

    return hcat(vertices...)'
end

function position_vectors_polygon(vertices, N, dist_type)
    V = length(vertices)
    all_x = Float64[]
    all_y = Float64[]
    for i in 1:V
        # Handle wrap-around at the last vertex
        p1 = vertices[i]
        p2 = vertices[i % V + 1]
        segment_points = interpolate_segment(p1, p2, N, dist_type)

        for i in 1:(length(segment_points) - 1)
            point = segment_points[i]
            push!(all_x, point[1])
            push!(all_y, point[2])
        end
    end
    return [all_x'; all_y']
end


function equidistant_points_on_polar_curve(x_function, y_function, num_points)

    function numerical_derivative(f, θ, h=1e-7)
        return (f(θ + h) - f(θ - h)) / (2h)
    end

    # Define the integrand for the arc length in polar coordinates
    integrand = θ -> sqrt(numerical_derivative(x_function, θ)^2 + numerical_derivative(y_function, θ)^2)

    function arc_length(θ)
        result, _ = quadgk(integrand, 0, θ)
        return result
    end

    # Equally spaced points along the polar curve in terms of arc length
    L, _ = quadgk(integrand, 0, 2π)  # Total length of the curve
    Δl = L / (num_points)
    Δθ = 2π / (num_points)

    theta_points = Float64[0.0]
    current_length = 0.0

   rootsFunc = (θ,curr_length) -> arc_length(θ) - (curr_length + Δl)


    for i in 1:num_points - 1
        θ = find_zero(θ->rootsFunc(θ,current_length), (theta_points[i], theta_points[i] + 2*Δθ))
        push!(theta_points, θ)
        current_length = arc_length(θ)
    end

    # Get equally spaced θ values along the polar curve
    θ_values = theta_points

    # Calculate corresponding (x, y) values
    x_values = x_function.(θ_values)
    y_values = y_function.(θ_values)

    return hcat(x_values, y_values)
end

"""
    u0SetUp(btype, R₀, N, dist_type)

Set up initial conditions for simulations based on the boundary type and distribution.

This function initializes the positions of particles or cells based on the specified boundary type and distribution.

# Arguments
- `btype`: Type of boundary (e.g., 'circle', 'triangle').
- `R₀`: Initial radius or characteristic length.
- `N`: Number of points or particles.
- `dist_type`: Type of distribution for the points.

# Returns
An array of initial positions.
"""
function u0SetUp(btype,R₀,N,dist_type,domain_type)
    # setting up initial conditions
    u0 = ElasticMatrix{Float64}(undef,2,N)

    if domain_type == "2D"
        if btype == "circle"
            R = R₀ # to produce identical areas
            θ = collect(NodeDistribution(0.0,2*π,N+1,dist_type)) 
            pop!(θ)
            @views u0 .= [X(R,θ)'; Y(R,θ)'];
        elseif btype == "triangle"
            #R = √((2*π*R₀^2)/sin(π/3))
            R = √((π*R₀^2)/(√(3)*cos(π/6)^2))
            # calc verticies
            vertices = polygon_vertices(3, R, -π/2)
            # calc number of nodes per segment 
            w = Int64(N/3) + 1
            @views u0 .= position_vectors_polygon(vertices, w, dist_type)
        elseif btype == "square"
            #R = √(π*(R₀^2)) # to produce identical areas
            R = (R₀√(2π))/2
            # calc verticies
            vertices = polygon_vertices(4, R, -π/4)
            # calc number of nodes per segment 
            w = Int64(N/4) + 1
            @views u0 .= position_vectors_polygon(vertices, w, dist_type)
        elseif btype == "hex"
            R = √((2/3√3)*π*(R₀^2)) # to produce identical areas
            # calc verticies
            vertices = polygon_vertices(6, R, 0)
            # calc number of nodes per segment 
            w = Int64(N/6) + 1
            @views u0 .= position_vectors_polygon(vertices, w, dist_type)
        elseif btype == "star"
            star_points = 5
            Rotation_Angle = pi/2
            rotation_angle = Rotation_Angle + pi/star_points
            vertices = StarVerticies(star_points, R₀, Rotation_Angle, rotation_angle)
            w = Int64(N/(2star_points)) + 1
            u0 .= position_vectors_polygon(vertices, w, dist_type)
        elseif btype == "cross"
            side_length = √((π*R₀^2)/5)
            offset = side_length/2
            vertices = CrossVertecies(side_length, offset)
            w = Int64(N/12) + 1
            @views u0 .= position_vectors_polygon(vertices, w, dist_type)
        elseif btype == "PerturbedCircle"
            Random.seed!(4) # nice ones: 333
            R_Pert = 2*R₀;
            x_range_loess = LinRange(0,2π,150)
            # generating random numbers
            ΔR = rand(150)
            # smoothing data with Loess
            loess_model = loess(x_range_loess,ΔR,span=0.3)
            θ_range = LinRange(0,2π,N)
            ΔR_loess = predict(loess_model, θ_range);
            # generating functions for equal distribution
            R = R₀ .+ ΔR_loess.*R_Pert
            θ = collect(NodeDistribution(0.0,2*π,N+1,dist_type)) 
            pop!(θ)
            @views u0 .= [X(R,θ)'; Y(R,θ)'];
            u0[:,end-5:end-1] = (u0[:,end-6:end-2] .+ u0[:,1]) ./ 2
            u0[:,end] = (u0[:,end-1] + u0[:,1]) / 2
        elseif btype == "CellCulture_P02"
            R = (2/3) * sqrt((8 * pi * (R₀^2) * (sin(deg2rad(36))^2)) / (sqrt(5*(5+2*sqrt(5))))) # to produce identical areas;
            h = R * cos(deg2rad(36));
            vertical_shift = (0, h)
            top_pentagon_vertices = polygon_vertices(5, R, π/2)
            bottom_pentagon_vertices = polygon_vertices(5, R, -π/2)
            vertices = [ top_pentagon_vertices[4] .+ vertical_shift, top_pentagon_vertices[5] .+ vertical_shift, top_pentagon_vertices[1] .+ vertical_shift, top_pentagon_vertices[2] .+ vertical_shift, top_pentagon_vertices[3] .+ vertical_shift,
                        bottom_pentagon_vertices[5] .- vertical_shift, bottom_pentagon_vertices[1] .- vertical_shift, bottom_pentagon_vertices[2] .- vertical_shift]
            w = Int64(N/(8)) + 1
            u0 .= position_vectors_polygon(vertices, w, dist_type)
        elseif btype == "VandenHeuvel_2023"
            R_img = [111.4398437,111.13292540811,110.558018488148,110.031317017362,109.548380134415,109.104766977972,108.696036686698,108.317748399258,107.965461254312,107.589591387666,107.294859447629,107.113860189539,106.82299376568,106.515206697927,106.171092909583,105.757472718751,105.335169354069,104.810978711499,104.262158336841,103.800047387094,103.276245009117,102.786551590443,102.291222952394,101.663100663949,101.076225621522,100.525271127923,99.8678357601493,99.075191901899,98.2524688643761,97.2436864012256,96.2620940734911,95.3689797926476,94.4208927691131,93.5911551417343,92.7475818629506,91.8041129570042,90.9274080087122,89.9727842667423,88.932299137363,87.8451718609023,86.7525826489298,85.6155235171774,84.5733348584536,83.4014210353829,82.3426257983279,81.4137241538084,80.4854517408829,79.600442595462,78.6831285154813,77.6992735211518,76.7015162188522,75.7063117826447,74.645386430636,73.4959568358492,72.410430451823,71.3960341835578,70.3896668751533,69.3861800615587,68.4902426833213,67.7262609626548,66.911694474329,66.0156568084702,65.0826222783384,64.2009021548825,63.2159776600115,62.1655332736809,61.0878486392854,60.0217535846467,59.0065787251315,58.0821022481891,57.2884936151888,56.6662550285924,56.1495491604789,55.8877575512315,55.8564881967475,55.8951628022683,56.0468581684912,56.4607695609945,57.1505191027715,57.9131612387337,58.7907869734484,59.8265681372757,60.9784483997073,62.1918933656603,63.4327250191094,64.69337993984,66.0050960530555,67.3598483649313,68.9554203105144,70.5388251521149,72.3084707040838,74.0980493226713,76.038523207695,77.8880501931432,79.7404753015705,81.5718288427254,83.4449575902097,85.244857465716,87.0297413169933,88.8463546384986,90.7315275009353,92.6006925691542,94.5505699878597,96.5766917415306,98.7199864727235,100.85613755025,102.911684698861,104.923813579984,107.011763963203,109.108899369898,110.693672366906,112.649269398314,114.995511481601,117.395691362097,120.304273488392,122.536698751117,124.376803981893,126.769016632257,128.539557392236,130.577662698114,131.981320564857,133.435268829961,135.067820035113,136.677915309266,138.378228783627,140.132553185211,142.793040214982,145.454812056181,148.233800806546,151.089163423997,153.477408456868,156.114751246763,159.035736536055,161.282727355307,163.350548396648,165.587023040557,167.778349376753,169.932169439657,172.025945172911,173.957194514561,175.741507125171,177.464916130205,179.212084768501,180.878182400425,182.269053482316,183.530484142689,184.704210749171,185.779719921919,186.658691245913,187.221572861429,187.340287335762,187.247036584589,187.007720319382,186.564102050042,185.943982436014,185.214766204204,184.489663879479,183.752424703818,182.879039239942,181.946194689849,181.191951863318,180.465235396646,179.769146583939,179.052630744294,178.264402115597,177.476353182275,176.873644043756,176.31092566497,175.648132815674,175.004922663644,174.291814689046,173.582063685018,172.792005211184,172.003136018723,171.186057845459,170.363385998641,169.573398106721,168.806372219188,167.974399758294,167.116608410353,166.242949942741,165.392259594668,164.567271002253,163.620505508806,162.621879113735,161.560762655815,160.42092604521,159.280677131148,158.176998778041,157.059569722829,155.927511962348,154.875645473768,153.794842994769,152.787661218058,151.933517884162,151.062704170058,150.18136924637,149.249639738529,148.265846926485,147.2901555882,146.253669977491,145.283349079058,144.39300367638,143.478016027565,142.457816588805,141.421916670894,140.518546884752,139.619705399027,138.896876137983,138.143385264794,137.3534882286,136.551634677557,135.797616672872,135.014452764789,134.160360339629,133.426332022687,132.625696475348,131.747681025653,130.755024833026,129.654318044004,128.538616912052,127.416811750579,126.354863164786,125.228495740338,124.05927022201,122.885377491737,121.69122651513,120.406598616078,119.071721315671,117.770832039039,116.464181574224,115.08130574873,113.695058853695,112.251668878677,110.840030958346,109.348420601541,107.880495246339,106.452211812232,105.066394315275,103.714980814464,102.445495126046,101.168973783606,99.8428590960953,98.6624969870973,97.5029714686967,96.4181738244076,95.2757617477299,94.243427444449,93.2645280687075,92.3048161359713,91.3214306453623,90.3781362625533,89.5750197750203,88.7313104046173,87.8641616458668,86.9878070909646,86.0172639185787,85.0309242300373,84.0322047452912,83.0166465898243,81.8699000883978,80.7855864569792,79.7480141212808,78.7501006185609,77.718517577941,76.7524269520809,75.8186944916924,74.9086746317886,73.9648486902811,73.0207573102727,72.0856076805104,71.2048058141613,70.2590603787801,69.2235415166968,68.0674677479408,66.9287446728471,65.755756445915,64.4871641888903,63.3515241852332,62.2336695574357,61.1449396612692,60.1264277246941,59.4252814921332,58.919727602316,58.7534601092257,58.821338512757,59.1888464922292,59.7628053031436,60.5295875951329,61.5123788649819,62.7083677561696,63.9791566350956,65.4116433623959,66.8771624074884,68.2495132226392,69.5827858085429,70.9964081887794,72.4956840277507,73.9940833360286,75.3665936872298,76.8634133653397,78.3571467640511,79.7845488320543,81.4904250143914,83.1525125629021,85.0443290708429,86.8621778550433,88.7153109934892,90.5916942045295,92.4383465467286,94.3003528100188,96.0834629179879,97.8590409069844,99.5782595742858,101.270616009186,103.156420936236,105.150295955302,106.872026597614,109.272648412408,111.182097641579,114.132148285118,116.125429052705,118.673505548862,121.215007748653,123.915947921685,126.087523868045,128.202383238625,129.853484351106,131.547342260087,133.104369280839,134.603578957452,136.200737378204,137.967653536786,139.829660030244,141.985705827418,144.161674140916,145.772933493566,148.632006963144,150.7655219824,153.224340042801,155.190507711782,156.596672195195,158.413868662934,159.835103295022,161.022100600784,161.921724386714,162.803609543211,162.988612990142,162.984459434621,162.652273484274,162.162095944517,161.395999695926,160.485797757588,159.490439936706,158.364565109046,157.198109623261,156.090655644701,155.069064764047,154.05660039799,153.025171518951,152.083830898786,151.214020673516,150.430530228051,149.643396023891,148.846359063799,148.055888822188,147.269627057281,146.481278535029,145.62220000747,144.660786849741,143.734386441884,142.777101731952,141.857619584539,140.941896748061,140.044156728184,139.216770818866,138.354619874654,137.483228233812,136.658077231548,135.870919654576,135.156708358375,134.421677672195,133.722812446713,132.938155988687,132.107401213825,131.304124738605,130.431042930619,129.503633937014,128.537458792671,127.54815811559,126.551448616953,125.669581809598,124.876317440096,124.041274721747,123.343644853753,122.623318737902,121.902536657601,121.132815324667,120.305871003253,119.4656030141,118.487173201316,117.524084731235,116.581793968047,115.669642535497,114.820423975006,114.055843520512,113.39760640595,112.867417865251,112.486983132351,112.278007441182,111.439843680699]
            θ_img = [0,0.00900876529041691,0.0270204491872647,0.045423279421577,0.0641321150014316,0.0823819305916906,0.101501831005953,0.119794022752755,0.139279259536846,0.159023179756261,0.177361810596894,0.195580679809826,0.208465218379987,0.226798848053885,0.24726176418276,0.267910422423339,0.288729528512397,0.309702944542456,0.321750554396642,0.346285507459378,0.367871670411153,0.380506377112364,0.40232109786044,0.415227335555231,0.441429043674087,0.459242351153995,0.477039665477051,0.49531123107642,0.509070888422383,0.519146114246522,0.547562235939997,0.561921562256815,0.573254410392972,0.593972681872105,0.608957619636508,0.621526624496621,0.639162741217661,0.652330901388891,0.665969237379109,0.688924388214861,0.703613782779452,0.718829999621624,0.743406056176078,0.768450633591043,0.785398163397448,0.80294022345485,0.828849058788979,0.847266009261828,0.866302262552678,0.885975080852296,0.914308935081439,0.935210619970345,0.956777289786955,0.979020156253393,1.00194846837353,1.02556917125762,1.04988654594668,1.07490183535883,1.10061286314815,1.12701365404525,1.15409406605549,1.18183944661588,1.21023032631679,1.23924216491895,1.26884516496464,1.29469930739946,1.32183362100217,1.35030256542822,1.38336721435122,1.41419444981288,1.42889927219073,1.46406065414576,1.4994888620096,1.53509721411557,1.57079632679489,1.60649543947422,1.64210379158018,1.67753199944403,1.71269338139906,1.74446752513647,1.77481430638744,1.80764508774181,1.83996381958059,1.86225312127276,1.88762081266311,1.91693229036425,1.94563704215049,1.96787577201987,1.98902065637412,2.01521553669599,2.0344439357957,2.04649154564988,2.0638472239997,2.08038859160382,2.08994244104141,2.11121582706548,2.12558855393279,2.13570072882567,2.14282584804053,2.16195349401094,2.16797298488757,2.1797539464314,2.19104581277771,2.20187571426622,2.21226904080412,2.22224955424568,2.2276554196435,2.23352583468788,2.236765564174,2.23983932498547,2.24968089171233,2.25527270552505,2.25774833097971,2.2655346029916,2.26869914938995,2.27402696783064,2.28119094367526,2.28287980208907,2.28527537490209,2.28888061439964,2.28962632641652,2.29035572088965,2.29245117765965,2.29504632754304,2.29626633507113,2.29685828448897,2.29743866747662,2.30310432923338,2.30410894097816,2.30507627128366,2.31035506431122,2.31549299871354,2.32049537751302,2.32112115965911,2.33011344890159,2.33055908167066,2.33510057293879,2.33545057817457,2.34384943834812,2.34404933839543,2.3521784475265,2.36017853285819,2.36805164199721,2.37579982104954,2.38342511171783,2.39092954863239,2.39831515690391,2.40558394988616,2.41273792713716,2.41672650997445,2.42831254653206,2.43933572208078,2.45025572482033,2.45849768683282,2.46685171136624,2.47531890952094,2.48390035778864,2.4948342273265,2.5067964391728,2.51561439514148,2.52134316760697,2.53359271606907,2.54275456080652,2.5551404363969,2.56451005278177,2.57702325992169,2.58956825844413,2.59922763048202,2.61188355194021,2.62173832416422,2.63170553400539,2.64178510244919,2.65469342177852,2.66495876166881,2.67533408446463,2.68581889761308,2.69641260183702,2.70711448760277,2.71792373171227,2.72883939404217,2.73986041445295,2.7509856098921,2.76221367171567,2.77354316325282,2.78282198331922,2.79440727235567,2.80609051029994,2.81587386605411,2.82781005805467,2.83603501571067,2.84828970736317,2.86063578206101,2.87134291014437,2.88393638126481,2.89661399046292,2.90784947272089,2.91930053516439,2.93237110105426,2.94551461603421,2.95749418232746,2.97085744211451,2.98211928414113,2.9947695304572,3.00857986671328,3.02160444351343,3.03560829560019,3.04900102378224,3.06318231485115,3.07693470716581,3.09091145351233,3.10511249446684,3.11953740689627,3.1342398449231,3.14899992552095,3.16381121891651,3.17916894293251,3.19457332949389,3.21018715704406,3.22600695926532,3.24202876563864,3.25824808903086,3.27465991559835,3.28460952855454,3.30278356802804,3.32000315294083,3.33898821343967,3.35669972672636,3.36839150164367,3.38657131671665,3.39891636856088,3.4198923125949,3.44136872368977,3.45498685653499,3.46051035437459,3.48339782557068,3.49503288394146,3.51009274657734,3.52557727754481,3.53271810055882,3.55786195144673,3.56578706149356,3.58302169726388,3.60972453746541,3.62796362343224,3.64175370733747,3.66073876783631,3.68021682365284,3.70019196893335,3.7116331172895,3.74164286699154,3.76311927808641,3.79392355497868,3.8075618909689,3.83203911064448,3.84651980572178,3.8701415696544,3.88569199861763,3.91872654232759,3.93539398053488,3.95262622550891,3.97044171237877,3.99705139962569,4.01604039962173,4.04357682244855,4.06374652871814,4.0845523865729,4.10600609400872,4.13545603536722,4.15808148418313,4.1813588755405,4.20529047599235,4.22987568536221,4.25511071200254,4.28098825381913,4.3074971940996,4.329104832639,4.357158118338,4.38083481850875,4.405712660983,4.43629196098925,4.46342627459197,4.49551536471828,4.52495986794101,4.55844431592634,4.59429695196559,4.62639468919044,4.66158527632646,4.69544145057828,4.72878095469269,4.76073836703659,4.79158805135414,4.82304620155858,4.84791669437019,4.86276040839182,4.8896143653344,4.91543419763815,4.93710014879933,4.96407895214287,4.98696426182431,5.00893478945439,5.03002433167422,5.05026716862496,5.06969766339278,5.08834993425599,5.10625759112105,5.12820076940931,5.13996976206931,5.15583731680442,5.17108613477665,5.1857450222196,5.19038967065931,5.19928821219595,5.20843265875086,5.21655494186384,5.22793798784366,5.23274428388468,5.24459848449856,5.2512084820206,5.26224216306127,5.268144701497,5.27844388917695,5.28840175082864,5.29091319480239,5.29336593840499,5.298169365368,5.30039158393225,5.30244919419716,5.30256549355118,5.30469263267914,5.30643989727462,5.30677447597178,5.30858812168216,5.31027268526032,5.31080787386322,5.31750018064446,5.31969420555094,5.32311494477389,5.32640801739263,5.33112331601179,5.3334071057242,5.34166077793879,5.34822740390843,5.35462905579848,5.36465081623049,5.37452147342323,5.38424253890948,5.39768615236417,5.41104880510615,5.42432649868856,5.43699258790734,5.44978910276963,5.45816032152821,5.47576427316916,5.48445460048799,5.49778714378213,5.50679590907255,5.52040771913745,5.52973973546661,5.54362656966325,5.55328564902785,5.56311461340549,5.57743052751856,5.59184837841012,5.60216829673246,5.61683871132782,5.63159812125672,5.64240158300207,5.65739351848372,5.66852235525742,5.67982882388961,5.69131425097559,5.7068100865884,5.71861591351148,5.73060095688887,5.742765806909,5.75511085875322,5.7676362997206,5.78383858549945,5.7966658702608,5.80966959745524,5.82284894435498,5.83620282189272,5.84972986305574,5.86342841175042,5.88027371632607,5.89133189935785,5.90834459182398,5.91989387204666,5.93441430359568,5.94651048779285,5.96143475278294,5.9765089877779,5.9894172687432,6.00488564817447,6.02049196186687,6.03623036204066,6.05209463998368,6.06807823404301,6.08253242985969,6.09732136135765,6.11388152902379,6.12924064272123,6.14494231670185,6.16204229302166,6.17830836844935,6.19568996637719,6.21250677312898,6.23013777741811,6.24780185105543,6.26533006204003,6.28318530717959]
            R_img_normalised = R_img ./ minimum(R_img)
            loess_model = loess(θ_img,R_img_normalised,span=0.1)
            θ_range = LinRange(0,2π,N);
            ΔR_loess = predict(loess_model, θ_range);
            R = R₀ .* ΔR_loess;
            θ = collect(NodeDistribution(0.0,2*π,N+1,dist_type)) 
            pop!(θ)
            @views u0 .= reverse!([X(R,θ)'; -Y(R,θ)'], dims=2);
        elseif btype == "Custom_Image"
            R_img, θ_img = KubaPhD.generate_r_θ_from_custom_image(readdir("./scripts/MPQC"; join=true)[1],500)
            #println(R_img, θ_img)
            R_img_normalised = R_img ./ minimum(R_img)
            loess_model = loess(θ_img,R_img_normalised,span=0.05)
            θ_range = LinRange(0,2π,N);
            ΔR_loess = predict(loess_model, θ_range);
            R = R₀ .* ΔR_loess;
            θ = collect(NodeDistribution(0.0,2*π,N+1,dist_type)) 
            pop!(θ)
            @views u0 .= [X(R,θ)'; -Y(R,θ)'];
        elseif btype == "Custom_Image1"
            R_img, θ_img = KubaPhD.generate_r_θ_from_custom_image(readdir("./scripts/MPQC"; join=true)[2],500)
            #println(R_img, θ_img)
            R_img_normalised = R_img ./ minimum(R_img)
            loess_model = loess(θ_img,R_img_normalised,span=0.05)
            θ_range = LinRange(0,2π,N);
            ΔR_loess = predict(loess_model, θ_range);
            R = R₀ .* ΔR_loess;
            θ = collect(NodeDistribution(0.0,2*π,N+1,dist_type)) 
            pop!(θ)
            @views u0 .= [X(R,θ)'; -Y(R,θ)'];
        end
    else
        if btype == "Line"
            xfunc = θ -> θ;
            yfunc = θ -> 0;
            @views u0 .= equidistant_points_on_polar_curve(xfunc, yfunc, N)';
        elseif btype == "SineWave"
            xfunc = θ -> θ;
            yfunc = θ -> 2 .+ 0.5 .* cos.(3 .* θ);
            #integrand(θ) = sqrt(numerical_derivative(xfunc, θ)^2 + numerical_derivative(yfunc, θ)^2)
            #rootsFunc(θ,curr_length,Δl) = arc_length(θ) - (curr_length + Δl)
            @views u0 .= equidistant_points_on_polar_curve(xfunc, yfunc, N)';
        elseif btype == "InvertedBellCurve"
            μ = 0.75 * 1000
            c = 0.2 * 1000^3
            xfunc = θ -> θ.*(1500/(2π));
            yfunc = θ -> -500 .* exp.((-((θ.*(1500/(2π))) .- μ).^6) ./ c.^2) .+ 500
            @views u0 .= equidistant_points_on_polar_curve(xfunc, yfunc, N)';
        end
    end
    return u0
end


# unchanged logic; slight cosmetic tidy
function color_mask(img::AbstractArray{<:Colorant}, which::Symbol)
    RGBimg = convert.(RGB{N0f8}, img)
    R = Float64.(channelview(RGBimg)[1, :, :])
    G = Float64.(channelview(RGBimg)[2, :, :])
    B = Float64.(channelview(RGBimg)[3, :, :])
    which === :red   && return (R .> 0.7) .& (G .< 0.3) .& (B .< 0.3)
    which === :green && return (G .> 0.7) .& (R .< 0.3) .& (B .< 0.3)
    error("which must be :red or :green")
end

function mask_boundary_xy(mask::AbstractMatrix{Bool})
    H, W = size(mask)
    z = transpose(Float64.(mask))        # (W, H) for Contour.jl
    x = collect(1:W)
    y = collect(1:H)

    cset = CTR.contours(x, y, z, [0.5])
    lvls = CTR.levels(cset)
    isempty(lvls) && return zeros(0,2)

    # pick the longest polyline from the single level
    best_line = nothing
    if size(lvls[1].lines, 1) > 1
        nmax = 0
        for ln in CTR.lines(lvls[1])
            n = length(CTR.coordinates(ln)[1])
            if n > nmax
                best_line = ln
                nmax = n
            end
        end
    else
        best_line = CTR.lines(lvls[1])[1]
    end
    best_line === nothing && return zeros(0,2)

    xs, ys = CTR.coordinates(best_line)
    return hcat(xs, ys)
end

"Resample a closed ring to exactly M points along arclength (open sequence, no duplicate endpoint)."
function resample_ring(xy::AbstractMatrix{<:Real}, M::Int)
    @assert size(xy,1) ≥ 3
    P = vcat(xy, xy[1, :]')                     # close
    d = sqrt.(sum(diff(P; dims=1).^2; dims=2))[:]
    s = cumsum(vcat(0.0, d)); total = s[end]
    st = range(0, total, length=M+1)[1:end-1]   # M samples, exclude duplicate end
    xt = similar(st); yt = similar(st)
    i = 1
    @inbounds for k in eachindex(st)
        t = st[k]
        while i < length(s)-1 && s[i+1] < t
            i += 1
        end
        α = (t - s[i]) / (s[i+1] - s[i] + eps())
        xt[k] = (1-α)*P[i,1] + α*P[i+1,1]
        yt[k] = (1-α)*P[i,2] + α*P[i+1,2]
    end
    return hcat(xt, yt)
end

# --- helpers for angles ---
function atan2pi(y, x)
    angle = atan(y, x)
    return angle >= 0 ? angle : angle + 2π
end

"Rotate a sequence so it starts at index k."
_rotate_start(A, k) = vcat(A[k:end, :], A[1:k-1, :])

"""
    generate_r_θ_from_custom_image(img_path::AbstractString, M::Int; closed::Bool=true)

Reads the red-mask boundary, centers it, resamples to M points, then returns a closed
polar sequence (`closed=true` ⇒ length M+1 with θ[1]=0 and θ[end]=2π).
Returns (r, θ, XY_centered). `r` and `θ` match in length (M or M+1).
"""
function generate_r_θ_from_custom_image(img_path::AbstractString, M::Int; closed::Bool=true)
    # centroid via polygon area formula (expects closed polygon; we close it locally)
    function calc_shape_centroid(xy)
        x = xy[:, 1]; y = xy[:, 2]
        P = vcat(xy, xy[1, :]')
        xP = P[:,1]; yP = P[:,2]
        cross = xP[1:end-1] .* yP[2:end] .- xP[2:end] .* yP[1:end-1]
        A = 0.5 * sum(cross)
        cx = (1 / (6A)) * sum((xP[1:end-1] .+ xP[2:end]) .* cross)
        cy = (1 / (6A)) * sum((yP[1:end-1] .+ yP[2:end]) .* cross)
        return [cx, cy]
    end

    img = load(img_path)
    r_xy = mask_boundary_xy(color_mask(img, :red))
    @assert size(r_xy,1) ≥ 3 "Could not extract a valid boundary from the red mask."

    C = calc_shape_centroid(r_xy)
    XY = r_xy .- C'                       # centered boundary (original density)
    XYs = resample_ring(XY, M)            # M points (open)

    # compute angles for resampled points in [0, 2π)
    θs = atan2pi.(XYs[:,2], XYs[:,1])

    # start the sequence at the smallest angle so θ[1] will become 0 after normalization
    k0 = argmin(θs)
    XYs = _rotate_start(XYs, k0)
    θs  = _rotate_start(hcat(θs, zeros(length(θs))), k0)[:,1]  # rotate θs consistently

    # normalize angles so first is exactly 0
    θs .-= θs[1]                 # now θs ∈ [0, 2π)
    r  = sqrt.(sum(XYs.^2; dims=2))[:]

    if closed
        # append the first point to close the loop
        XYc = vcat(XYs, XYs[1, :]')
        rc  = vcat(r, r[1])
        θc  = vcat(θs, 2π)        # force exact 2π at the last sample
        return rc, θc, XY
    else
        # open sequence spanning [0, 2π) (no explicit 2π point)
        return r, θs, XY
    end
end


"""
    NodeDistribution(start, stop, length, type)

Generate a distribution of nodes between `start` and `stop` with a specified distribution type.

# Arguments
- `start`: Starting value of the range.
- `stop`: Ending value of the range.
- `length`: Number of points in the distribution.
- `type`: Type of distribution (e.g., 'Linear', 'exp', 'sine').

# Returns
A range or array of distributed values.
"""
function NodeDistribution(start,stop,length,type)
    if type == "Linear"
        return LinRange(start, stop, length)
    else
        return nonLinearRange(start, stop, length, type)
    end
end

function nonLinearRange(start, stop, length, dist_type)
    linear_range = LinRange(0, 1, length)

    # Applying different distribution types
    if dist_type == "exp"
        # Exponential scaling
        return start .+ (exp.(linear_range .* log(1 + stop - start)) .- 1)
    elseif dist_type == "sine"
        # Sinusoidal scaling
        return start .+ (sin.((π/2) .* linear_range) .* (stop - start))
    elseif dist_type == "cosine"
        # Cosine scaling
        return start .+ ((1 .- cos.((π/2) .* linear_range)) .* (stop - start))
    elseif dist_type == "quad"
        # Quadratic scaling
        return start .+ (linear_range .^ 2 .* (stop - start))
    elseif dist_type == "cubic"
        # Cubic scaling
        return start .+ (linear_range .^ 3 .* (stop - start))
    elseif dist_type == "sigmoid"
        # Sigmoid scaling
        linear_range = LinRange(-1, 1, length)
        k = 5; # slope steepness
        sigmoid_range = 1 ./ (1 .+ exp.(-k.*(linear_range)))
        return start .+ (sigmoid_range .* (stop - start))
    elseif dist_type == "2sigmoid"
        # Piecewise sigmoid scaling
        k1 = 10;  k2 = 10;
        x01 = 0.5;  x02 = 0.5;
        piecewise_sigmoid = [x < 0.5 ? 0.5 * (1 / (1 + exp(-k1 * (2x - x01)))) : 0.5 + 0.5 * (1 / (1 + exp(-k2 * (2x - 1 - x02)))) for x in linear_range]
        return start .+ (piecewise_sigmoid * (stop - start))
    elseif dist_type == "4sigmoid"
        # Parameters for the sigmoid functions
        k = 20
        # Adjust midpoints for the full sigmoid in the first and last segments
        x0 = [0.125, 0.375, 0.625, 0.875]

        # Piecewise sigmoid scaling with 4 segments
        piecewise_sigmoid = [if x < 0.25
                                (1 / (1 + exp(-k * (4x - x0[1])))) * 0.25
                             elseif x < 0.5
                                0.25 + (1 / (1 + exp(-k * (4x - 1 - x0[2])))) * 0.25
                             elseif x < 0.75
                                0.5 + (1 / (1 + exp(-k * (4x - 2 - x0[3])))) * 0.25
                             else
                                0.75 + (1 / (1 + exp(-k * (4x - 3 - x0[4])))) * 0.25
                             end for x in linear_range]

        return start .+ (piecewise_sigmoid .* (stop - start))
    else
        error("Unsupported distribution type")
    end
end