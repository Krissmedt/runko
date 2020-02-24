
# all headers
# python3 tools/run-clang-tidy.py -p=build/ \
#	-header-filter=.* \
#	-checks=-*,modernize-* 

# pick only runko headers
python3 tools/run-clang-tidy.py -p=build/ \
	-header-filter=/home/natj/runko \
	-fix \
	-checks=-*,\
modernize-avoid-bind,\
modernize-concat-nested-namespaces,\
modernize-deprecated-headers,\
modernize-deprecated-ios-base-aliases,\
modernize-make-shared,\
modernize-make-unique,\
modernize-pass-by-value,\
modernize-raw-string-literal,\
modernize-redundant-void-arg,\
modernize-replace-auto-ptr,\
modernize-replace-random-shuffle
#modernize-return-braced-init-list,\
#modernize-shrink-to-fit,\
#modernize-unary-static-assert,\
#modernize-use-auto,\
#modernize-use-default-member-init,\
#modernize-use-bool-literals,\
#modernize-use-emplace,\
#modernize-use-equals-default,\
#modernize-use-equals-delete,\
#modernize-use-nodiscard,\
#modernize-use-noexcept

# -checks=-*,modernize-*,-modernize-use-trailing-return-type,-modernize-avoid-c-arrays

#checks not to be run
# modernize-loop-convert


# modernize-pass-by-value
# modernize-use-equals-default
# modernize-use-default-member-init
# modernize-use-using

# modernize-deprecated-headers
	#-checks=-*,modernize-*,-modernize-use-trailing-return-type 


	#-checks=-*,modernize-avoid-bind

