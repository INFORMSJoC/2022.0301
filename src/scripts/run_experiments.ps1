
function run_experiments {

	Param(
		[Parameter(Mandatory)] 
		[string] $n,
		[Parameter(Mandatory)]
		[string] $m)
	
	Get-ChildItem .\instances\ | Where-Object { $_.Name -match "wt($n)" } | ForEach-Object { Get-ChildItem "$_" | Where-Object { $_.Name -match ".*[16]\.dat" } | ForEach-Object -Parallel {.\build\PM.exe  -j .\settings\default_settings.json "$_" "$using:m"}  }  
}
