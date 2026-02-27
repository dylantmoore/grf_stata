# Makefile for grf_plugin (C++ wrapper around grf library)
#
# Targets:
#   all            - Build for local platform (macOS ARM64)
#   macosx         - Build for macOS Apple Silicon (ARM64)
#   linux          - Cross-compile for Linux x86_64
#   windows        - Cross-compile for Windows x86_64
#   all-platforms  - Build all three targets
#   clean          - Remove all built plugins

VENDOR_GRF = vendor/grf/core
GRF_CORE   = $(VENDOR_GRF)/src
GRF_3RDPTY = $(VENDOR_GRF)/third_party

INCLUDES = -I$(GRF_CORE) -I$(GRF_3RDPTY) -I. -Wno-deprecated-declarations

# All grf C++ source files
GRF_SRCS = \
    $(GRF_CORE)/commons/Data.cpp \
    $(GRF_CORE)/commons/utility.cpp \
    $(GRF_CORE)/forest/Forest.cpp \
    $(GRF_CORE)/forest/ForestOptions.cpp \
    $(GRF_CORE)/forest/ForestPredictor.cpp \
    $(GRF_CORE)/forest/ForestPredictors.cpp \
    $(GRF_CORE)/forest/ForestTrainer.cpp \
    $(GRF_CORE)/forest/ForestTrainers.cpp \
    $(GRF_CORE)/prediction/Prediction.cpp \
    $(GRF_CORE)/prediction/PredictionValues.cpp \
    $(GRF_CORE)/prediction/RegressionPredictionStrategy.cpp \
    $(GRF_CORE)/prediction/InstrumentalPredictionStrategy.cpp \
    $(GRF_CORE)/prediction/MultiCausalPredictionStrategy.cpp \
    $(GRF_CORE)/prediction/MultiRegressionPredictionStrategy.cpp \
    $(GRF_CORE)/prediction/ProbabilityPredictionStrategy.cpp \
    $(GRF_CORE)/prediction/QuantilePredictionStrategy.cpp \
    $(GRF_CORE)/prediction/SurvivalPredictionStrategy.cpp \
    $(GRF_CORE)/prediction/CausalSurvivalPredictionStrategy.cpp \
    $(GRF_CORE)/prediction/LocalLinearPredictionStrategy.cpp \
    $(GRF_CORE)/prediction/LLCausalPredictionStrategy.cpp \
    $(GRF_CORE)/prediction/ObjectiveBayesDebiaser.cpp \
    $(GRF_CORE)/prediction/collector/DefaultPredictionCollector.cpp \
    $(GRF_CORE)/prediction/collector/OptimizedPredictionCollector.cpp \
    $(GRF_CORE)/prediction/collector/SampleWeightComputer.cpp \
    $(GRF_CORE)/prediction/collector/TreeTraverser.cpp \
    $(GRF_CORE)/relabeling/NoopRelabelingStrategy.cpp \
    $(GRF_CORE)/relabeling/MultiNoopRelabelingStrategy.cpp \
    $(GRF_CORE)/relabeling/InstrumentalRelabelingStrategy.cpp \
    $(GRF_CORE)/relabeling/MultiCausalRelabelingStrategy.cpp \
    $(GRF_CORE)/relabeling/QuantileRelabelingStrategy.cpp \
    $(GRF_CORE)/relabeling/LLRegressionRelabelingStrategy.cpp \
    $(GRF_CORE)/relabeling/CausalSurvivalRelabelingStrategy.cpp \
    $(GRF_CORE)/sampling/RandomSampler.cpp \
    $(GRF_CORE)/sampling/SamplingOptions.cpp \
    $(GRF_CORE)/splitting/RegressionSplittingRule.cpp \
    $(GRF_CORE)/splitting/InstrumentalSplittingRule.cpp \
    $(GRF_CORE)/splitting/MultiCausalSplittingRule.cpp \
    $(GRF_CORE)/splitting/MultiRegressionSplittingRule.cpp \
    $(GRF_CORE)/splitting/ProbabilitySplittingRule.cpp \
    $(GRF_CORE)/splitting/SurvivalSplittingRule.cpp \
    $(GRF_CORE)/splitting/AcceleratedSurvivalSplittingRule.cpp \
    $(GRF_CORE)/splitting/CausalSurvivalSplittingRule.cpp \
    $(GRF_CORE)/splitting/factory/RegressionSplittingRuleFactory.cpp \
    $(GRF_CORE)/splitting/factory/InstrumentalSplittingRuleFactory.cpp \
    $(GRF_CORE)/splitting/factory/MultiCausalSplittingRuleFactory.cpp \
    $(GRF_CORE)/splitting/factory/MultiRegressionSplittingRuleFactory.cpp \
    $(GRF_CORE)/splitting/factory/ProbabilitySplittingRuleFactory.cpp \
    $(GRF_CORE)/splitting/factory/SurvivalSplittingRuleFactory.cpp \
    $(GRF_CORE)/splitting/factory/CausalSurvivalSplittingRuleFactory.cpp \
    $(GRF_CORE)/tree/Tree.cpp \
    $(GRF_CORE)/tree/TreeOptions.cpp \
    $(GRF_CORE)/tree/TreeTrainer.cpp \
    $(GRF_CORE)/analysis/SplitFrequencyComputer.cpp \
    $(GRF_CORE)/commons/ProgressBar.cpp \
    $(GRF_CORE)/RuntimeContext.cpp

# Plugin source
PLUGIN_SRC = grf_plugin.cpp

# stplugin.c — compiled as C separately per platform
STPLUGIN_SRC = stplugin.c

# Output filenames
TARGET_DARWIN_ARM64  = grf_plugin_macosx.plugin
TARGET_LINUX         = grf_plugin_unix.plugin
TARGET_WINDOWS       = grf_plugin_windows.plugin

# ── Darwin arm64 (macOS Apple Silicon) ─────────────────────────────
DARWIN_ARM64_CXX    = g++
DARWIN_ARM64_CC     = gcc
DARWIN_ARM64_CXXFLAGS = -std=c++17 -O3 -fPIC -DSYSTEM=APPLEMAC -arch arm64 $(INCLUDES)
DARWIN_ARM64_CFLAGS   = -O3 -fPIC -DSYSTEM=APPLEMAC -arch arm64 -I.
DARWIN_ARM64_LDFLAGS  = -bundle -lpthread

# ── Windows (x86_64, cross-compiled with mingw-w64) ───────────────
WIN_CXX    = x86_64-w64-mingw32-g++
WIN_CC     = x86_64-w64-mingw32-gcc
WIN_CXXFLAGS = -std=c++17 -O3 -DSYSTEM=STWIN32 $(INCLUDES)
WIN_CFLAGS   = -O3 -DSYSTEM=STWIN32 -I.
WIN_LDFLAGS  = -shared -static-libstdc++ -static-libgcc -lpthread

# ── Linux (x86_64) ────────────────────────────────────────────────
LINUX_CXX    = g++
LINUX_CC     = gcc
LINUX_CXXFLAGS = -std=c++17 -O3 -fPIC -DSYSTEM=OPUNIX $(INCLUDES)
LINUX_CFLAGS   = -O3 -fPIC -DSYSTEM=OPUNIX -I.
LINUX_LDFLAGS  = -shared -static-libstdc++ -static-libgcc -lpthread

# ── Phony targets ─────────────────────────────────────────────────
.PHONY: all macosx linux windows all-platforms clean

# Default: build for local platform only
all: macosx

macosx: $(TARGET_DARWIN_ARM64)
linux: $(TARGET_LINUX)
windows: $(TARGET_WINDOWS)
all-platforms: macosx linux windows

# ── Build rules ───────────────────────────────────────────────────

$(TARGET_DARWIN_ARM64): $(PLUGIN_SRC) $(GRF_SRCS) $(STPLUGIN_SRC)
	$(DARWIN_ARM64_CC) $(DARWIN_ARM64_CFLAGS) -c $(STPLUGIN_SRC) -o stplugin.darwin-arm64.o
	$(DARWIN_ARM64_CXX) $(DARWIN_ARM64_CXXFLAGS) $(DARWIN_ARM64_LDFLAGS) -o $@ $(PLUGIN_SRC) $(GRF_SRCS) stplugin.darwin-arm64.o
	rm -f stplugin.darwin-arm64.o

$(TARGET_LINUX): $(PLUGIN_SRC) $(GRF_SRCS) $(STPLUGIN_SRC)
	$(LINUX_CC) $(LINUX_CFLAGS) -c $(STPLUGIN_SRC) -o stplugin.linux.o
	$(LINUX_CXX) $(LINUX_CXXFLAGS) $(LINUX_LDFLAGS) -o $@ $(PLUGIN_SRC) $(GRF_SRCS) stplugin.linux.o
	rm -f stplugin.linux.o

$(TARGET_WINDOWS): $(PLUGIN_SRC) $(GRF_SRCS) $(STPLUGIN_SRC)
	$(WIN_CC) $(WIN_CFLAGS) -c $(STPLUGIN_SRC) -o stplugin.windows.o
	$(WIN_CXX) $(WIN_CXXFLAGS) $(WIN_LDFLAGS) -o $@ $(PLUGIN_SRC) $(GRF_SRCS) stplugin.windows.o
	rm -f stplugin.windows.o

clean:
	rm -f $(TARGET_DARWIN_ARM64) $(TARGET_LINUX) $(TARGET_WINDOWS)
	rm -f stplugin.*.o
