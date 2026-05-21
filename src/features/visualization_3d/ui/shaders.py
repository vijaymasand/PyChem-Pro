"""
GLSL Shaders for high-performance molecular visualization.
Optimized for high-end GPUs like the RTX 3060.
"""

# ─── LINE SHADER (Bonds/Sticks) ───────────────────────────────────────────────
LINE_VERTEX_SHADER = """
#version 330 core
layout(location = 0) in vec3 aPos;
layout(location = 1) in vec3 aColor;

uniform mat4 view;
uniform mat4 projection;

out vec3 vColor;

void main() {
    vColor = aColor;
    vec4 clipPos = projection * view * vec4(aPos, 1.0);
    // Pull bonds forward significantly so they aren't hidden inside spheres
    clipPos.z -= 0.025 * clipPos.w; 
    gl_Position = clipPos;
}
"""

LINE_FRAGMENT_SHADER = """
#version 330 core
out vec4 FragColor;
in vec3 vColor;

void main() {
    FragColor = vec4(vColor, 1.0);
}
"""

# ─── MESH SHADER (Cartoon/Ribbon) ─────────────────────────────────────────────
MESH_VERTEX_SHADER = """
#version 330 core
layout(location = 0) in vec3 aPos;
layout(location = 1) in vec3 aColor;
layout(location = 2) in vec3 aNormal;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

out vec3 FragPos;
out vec3 Normal;
out vec3 Color;

void main() {
    FragPos = vec3(model * vec4(aPos, 1.0));
    Normal = mat3(transpose(inverse(model))) * aNormal;
    Color = aColor;
    gl_Position = projection * view * vec4(FragPos, 1.0);
}
"""

MESH_FRAGMENT_SHADER = """
#version 330 core
out vec4 FragColor;

in vec3 FragPos;
in vec3 Normal;
in vec3 Color;

uniform vec3 lightPos;
uniform vec3 viewPos;
uniform vec3 lightColor;

void main() {
    // Ambient
    float ambientStrength = 0.35;
    vec3 ambient = ambientStrength * lightColor;
    
    // Diffuse
    vec3 norm = normalize(Normal);
    vec3 lightDir = normalize(lightPos - FragPos);
    float diff = max(dot(norm, lightDir), 0.0);
    vec3 diffuse = diff * lightColor;
    
    // Specular
    float specularStrength = 0.5;
    vec3 viewDir = normalize(viewPos - FragPos);
    vec3 reflectDir = reflect(-lightDir, norm);
    float spec = pow(max(dot(viewDir, reflectDir), 0.0), 32);
    vec3 specular = specularStrength * spec * lightColor;
    
    vec3 result = (ambient + diffuse + specular) * Color;
    FragColor = vec4(result, 1.0);
}
"""

# ─── SPHERE IMPOSTOR SHADER (Atoms) ──────────────────────────────────────────
# Renders a screen-aligned quad as a perfect shaded sphere.
# Fast, pixel-perfect, and handles depth correctly.
SPHERE_VERTEX_SHADER = """
#version 330 core
layout(location = 0) in vec3 aPos;     // Center of sphere
layout(location = 1) in vec3 aColor;   // Color
layout(location = 2) in float aRadius; // Radius

uniform mat4 view;
uniform mat4 projection;
uniform float zoom;

out vec3 vCenter;
out vec3 vColor;
out float vRadius;
out vec2 vTexCoord;

void main() {
    vCenter = aPos;
    vColor = aColor;
    vRadius = aRadius;
    
    // Calculate billboard offsets
    vec3 cameraRight = vec3(view[0][0], view[1][0], view[2][0]);
    vec3 cameraUp = vec3(view[0][1], view[1][1], view[2][1]);
    
    // Quad vertices around center (6 per quad)
    vec2 offsets[6] = vec2[](
        vec2(-1,-1), vec2(1,-1), vec2(1,1),
        vec2(-1,-1), vec2(1,1), vec2(-1,1)
    );
    vTexCoord = offsets[gl_VertexID % 6];
    
    vec3 pos = aPos + (cameraRight * vTexCoord.x + cameraUp * vTexCoord.y) * vRadius;
    gl_Position = projection * view * vec4(pos, 1.0);
}
"""

SPHERE_FRAGMENT_SHADER = """
#version 330 core
out vec4 FragColor;

in vec3 vCenter;
in vec3 vColor;
in float vRadius;
in vec2 vTexCoord;

uniform mat4 projection;
uniform mat4 view;
uniform vec3 lightPos;

void main() {
    float r2 = dot(vTexCoord, vTexCoord);
    if (r2 > 1.0) discard; // Hard cut
    
    // Soften edge slightly to help aliasing
    float edge = 1.0 - smoothstep(0.96, 1.0, r2);
    
    // Calculate 3D normal from 2D coordinates
    float z = sqrt(1.0 - r2);
    vec3 normal = vec3(vTexCoord, z);
    
    // Shading
    vec3 lightDir = normalize(vec3(0.3, 0.5, 1.0)); 
    float diff = max(dot(normal, lightDir), 0.0);
    float spec = pow(max(dot(normal, normalize(lightDir + vec3(0,0,1))), 0.0), 40.0);
    
    vec3 ambient = 0.4 * vColor;
    vec3 diffuse = 0.6 * diff * vColor;
    vec3 specular = 0.3 * spec * vec3(1.0);
    
    FragColor = vec4(ambient + diffuse + specular, edge);
    
    // Update depth buffer
    vec4 viewPos = view * vec4(vCenter, 1.0);
    viewPos.z += z * vRadius; 
    vec4 clipPos = projection * viewPos;
    gl_FragDepth = (clipPos.z / clipPos.w) * 0.5 + 0.5;
}
"""
