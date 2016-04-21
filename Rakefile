@filename = "CutTheRoots"

task :default do
  `g++ -std=c++11 -W -Wall -Wno-sign-compare -O2 -o #{@filename} #{@filename}.cpp`
end

task :run do
  system("g++ -W -Wall -Wno-sign-compare -O2 -o #{@filename} #{@filename}.cpp")
  system("java -jar ./visualizer.jar -size 8 -vis -save result.png -seed 1 -exec './#{@filename}'")
  #system("java -jar ./#{@filename}Vis.jar -side 12 -seed 105 -exec './#{@filename}'")
end

task :windows do
  system("g++ -std=c++11 -W -Wall -Wno-sign-compare -O2 -o #{@filename} #{@filename}.cpp")
  system("java -jar ./visualizer.jar -novis -save result.png -seed 10 -exec ./#{@filename}.exe")
end

task :one do
  system("g++ -std=c++11 -W -Wall -Wno-sign-compare -O2 -o #{@filename} #{@filename}.cpp")
  system("time java -jar visualizer.jar -seed 2987 -save result.png -novis -exec './#{@filename}'")
end

task :debug do
  system("g++ -std=c++11 -W -Wall -g -fsanitize=address -fno-omit-frame-pointer -Wno-sign-compare -O2 -o #{@filename} #{@filename}.cpp")
  system("time java -jar visualizer.jar -seed 3 -save result.png -novis -exec './#{@filename}'")
end

task :cover do
  system("g++ -W -Wall -Wno-sign-compare -o #{@filename} --coverage #{@filename}.cpp")
  system("time java -jar visualizer.jar -seed 10 -save result.png -novis -exec './#{@filename}'")
end

task :clean do
  system("rm *.gcda")
  system("rm *.gcov")
  system("rm *.gcno")
end

task :two do
  system("g++ -std=c++11 -W -Wall -Wno-sign-compare -O2 -o #{@filename} #{@filename}.cpp")
  system("time java -jar visualizer.jar -seed 10 -novis -exec './#{@filename}'")
end

task :novis do
  system('rm result.txt')
  system("g++ -W -Wall -Wno-sign-compare -O2 -o #{@filename} #{@filename}.cpp")
  1001.upto(1100) do |num|
    p num
    system("time java -jar ./visualizer.jar -seed #{num} -novis -exec './#{@filename}' >> result.txt")
  end
  system('ruby analysis.rb 100')
end

task :final do
  system('rm result.txt')
  system("g++ -std=c++11 -W -Wall -Wno-sign-compare -O2 -o #{@filename} #{@filename}.cpp")
  2001.upto(3000) do |num|
    p num
    system("time java -jar ./visualizer.jar -seed #{num} -novis -exec './#{@filename}' >> result.txt")
  end
  system('ruby analysis.rb 1000')
end

task :sample do
  system('rm result.txt')
  system("g++ -std=c++11 -W -Wall -Wno-sign-compare -O2 -o #{@filename} #{@filename}.cpp")
  1.upto(10) do |num|
    p num
    system("time java -jar ./visualizer.jar -seed #{num} -novis -exec './#{@filename}' >> result.txt")
  end
  system('ruby analysis.rb 10')
end

task :select do
  system('rm result.txt')
  system("g++ -W -Wall -Wno-sign-compare -O2 -o #{@filename} #{@filename}.cpp")
  array = [2163, 2175, 2184, 2227, 2249, 2272, 2307, 2320, 2381, 2391, 2406, 2444, 2493, 2495, 2504, 2509, 2510, 2515, 2521, 2552, 2553,
    2566, 2573, 2590, 2597, 2599, 2601, 2603, 2615, 2617, 2632, 2647, 2671, 2675, 2676, 2685, 2712, 2718, 2751, 2752, 2776, 2787, 2804,
    2808, 2817, 2819, 2839, 2846, 2850, 2856, 2862, 2863, 2877, 2893, 2918, 2927, 2934, 2950, 2953, 2997, ]
  array.take(100).each do |num|
    p num
    system("time java -jar ./visualizer.jar -seed #{num} -novis -exec './#{@filename}' >> result.txt")
  end
  system("ruby analysis.rb #{array.size}")
end

task :test do
  system("g++ -o #{@filename} #{@filename}.cpp")
  system("./#{@filename} < test_case.txt")
end

