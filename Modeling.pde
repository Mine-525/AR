


void draw_cylinder(float top_radius, float bottom_radius, float height){
    int sides = 20;
    float angle = 0;
    float angle_increment = radians(360 / 20);

    // side
    beginShape(QUAD_STRIP);
        for(int i = 0; i < sides + 1; i++){
            vertex(top_radius * cos(angle), height, top_radius * sin(angle));
            vertex(bottom_radius * cos(angle), 0, bottom_radius * sin(angle));
            angle += angle_increment;
        }
    endShape();


    angle = 0;
    beginShape(TRIANGLE_FAN);
        // center point
        vertex(0, height, 0);
        for (int i = 0; i < sides + 1; i++){
            vertex(bottom_radius * cos(angle), 0, bottom_radius * sin(angle));
            angle += angle_increment;
        }
    endShape();
    

    if(top_radius != 0){
        angle = 0;
        beginShape(TRIANGLE_FAN);
            // center point
            vertex(0, 0, 0);
            for (int i = 0; i < sides + 1; i++){
                vertex(top_radius * cos(angle), height, top_radius * sin(angle));
                angle += angle_increment;
            }
        endShape();
    }
}

void snowman(float scale){
    float bottom_radius = 2 * scale;
    float mid_radius = bottom_radius * 0.75;
    float top_radius = bottom_radius * 0.5;

    // origin is the center of the bottom sphere
    // body
    // specular(100);
    noStroke();
    fill(255);
    sphere(bottom_radius);

    translate(0, -bottom_radius * 1.25, 0);

    sphere(mid_radius);

    translate(0, -bottom_radius, 0);

    // face
    sphere(top_radius);

    // eye
    float theta = 30;
    float phi = 20;
    float eye_radius = 0.1;
    float r = top_radius - eye_radius;
    // specular(0, 0, 0);
    // shininess(0);
    fill(0);

    pushMatrix();
        translate(-r * cos(radians(phi)) * sin(radians(theta)), 
                -r * sin(radians(phi)), 
                r * cos(radians(phi) * cos(radians(theta))));
        ambient(0);
        sphere(0.1);
    popMatrix();

    pushMatrix();
        translate(r * cos(radians(phi)) * sin(radians(theta)), 
                -r * sin(radians(phi)), 
                r * cos(radians(phi) * cos(radians(theta))));
        sphere(0.1);
    popMatrix();

    // nose
    fill(255, 120, 0);
    theta = 0;
    phi = 10;
    r -= 0.05;
    pushMatrix();
        translate(r * cos(radians(phi)) * sin(radians(theta)), 
                -r * sin(radians(phi)), 
                r * cos(radians(phi) * cos(radians(theta))));
        
        rotateX(radians((90 + phi)));

        ambient(255, 120, 0);
        draw_cylinder(0, 0.2, 0.5);
    popMatrix();





}
