submodule (maptran) no_sphere
implicit none

contains

module procedure lookAtSpheroid
  error stop 'this compiler does not have working NaN'
end procedure lookAtSpheroid

end submodule no_sphere