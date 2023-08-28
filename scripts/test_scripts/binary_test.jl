function minimize_binary_search(list::Vector{Float64}, target::Float64)
    left = 1
    right = length(list)
    while left < right
        mid = div(left + right,  2)
        if list[mid] < target
            # list[mid + 1] < target ? (left = mid + 1) : return mid, mid + 1
            left = mid + 1
        elseif list[mid] > target
            # list[mid - 1] > target ? (right = mid - 1) : return mid - 1, mid
            right = mid - 1
        else
            return mid, mid
        end 
    end
    return left, left + 1
end

function main()
    a = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]
    for i in 1.0:0.1:9.0
        @show minimize_binary_search(a, i)
    end
end

main()