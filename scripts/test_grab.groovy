
@GrabResolver(name='central', root='https://repo1.maven.org/maven2/')
@Grab(group='com.github.sharispe', module='slib-sml', version='0.9.1')
import slib.sml.sm.core.metrics.utils.SMConstants
println "SML SMConstants loaded: " + SMConstants.FLAG_SIM_GROUP_BMA
